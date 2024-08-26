use crate::grn::grn_protein;

use super::grn_protein::GrnProtein;
// use jaded::{Content, Parser, Result};
use std::{fs::File, io, cmp};
use serde::Deserialize;

#[derive(Deserialize, Debug)]
#[allow(non_snake_case)]
struct DataReaded {
    proteins : String,
    beta : f64,
    delta : f64,
    maxEnhance : f64,
    maxInhibit : f64,
}

#[derive(Debug)]
pub struct GrnModel{ 
    beta: f64,
    delta: f64,
    current_step: i32,
    max_enhance: f64,
    max_inhibit: f64,
    enhance_matching: Vec<Vec<f64>>,
    inhibit_matching: Vec<Vec<f64>>,
    pub proteins: Vec<GrnProtein>,
}


impl GrnModel{

    // constructor with parameters
    pub fn new(p : Vec<GrnProtein>, beta : f64, delta : f64) -> GrnModel
    {
        GrnModel {
            // initialisation of proteins
            beta : beta,
            delta : delta,
            proteins : p,
            current_step : 0,
            max_enhance : 0.,
            max_inhibit : 0.,
            enhance_matching : Vec::new(),
            inhibit_matching : Vec::new()
        }
    }

    // constructor without parameters
    pub fn new_default() -> GrnModel
    {
        GrnModel {
            // initialisation of proteins
            beta : 0.0,
            delta : 0.0,
            proteins : Vec::new(),
            current_step : 0,
            max_enhance : 0.,
            max_inhibit : 0.,
            enhance_matching : Vec::new(),
            inhibit_matching : Vec::new()
        }
    }

    pub fn load_from_file(&mut self, filename : &str) {
    
        // let sample = File::open(filename).expect("File missing");
        // let mut parser = Parser::new(sample);
        // let obj: std::result::Result<jaded::Content, jaded::JavaError> = parser.unwrap().read();
        // println!("Readed object {:?}", obj.unwrap().value());
        // // TODO implement the reading strategy to communicate between java and rust
    
        let sample = File::open(filename).expect("File missing");
        let reader = io::BufReader::new(sample);
        
        // Deserialize the JSON data into the struct
        let data: DataReaded = serde_json::from_reader(reader).expect("custom");
        //println!("readed data is ! {:?}", data);
        
        self.beta = data.beta;
        self.delta = data.delta;
        self.max_enhance = data.maxEnhance;
        self.max_inhibit = data.maxInhibit;
        
        // convert proteins from_string;
        let proteins : Vec<&str>= data.proteins.split("][").collect();

        for p in proteins
        {
            let temp_prot = grn_protein::GrnProtein::from_string(p.to_string());
            self.proteins.push(temp_prot);
        }

        self.update_signature()
        
    }

    fn update_signature(&mut self)
    { // function that compute the max and matching of  inhibiting and enhancing
        let prot_size = self.proteins.len();

        // compute max values
        let mut enh_match = vec![vec![0.0; prot_size]; prot_size];
        let mut inh_match = vec![vec![0.0; prot_size]; prot_size];
        let mut max_enh: f64 = 0.;
        let mut max_inh: f64 = 0.;


        for j in 0..prot_size
        {
            
            for k in 0..prot_size
            {
				enh_match[j][k] = (grn_protein::IDSIZE-(self.proteins[j].enhancer-self.proteins[k].id).abs() ) as f64;
				max_enh= max_enh.max(enh_match[j][k]);
				inh_match[j][k] = (grn_protein::IDSIZE-(self.proteins[j].inhibiter-self.proteins[k].id).abs() ) as f64;
				max_inh= max_inh.max(inh_match[j][k]);                
            }
        }

        for j in 0..prot_size
        {
            
            for k in 0..prot_size
            {
				enh_match[j][k] = (self.beta * enh_match[j][k] - max_enh).exp();			    
                inh_match[j][k] = (self.beta * inh_match[j][k] - max_inh).exp(); 
			}
        }

        self.max_enhance = max_enh;
        self.max_inhibit = max_inh;
        self.enhance_matching = enh_match;
        self.inhibit_matching = inh_match;
        
    }   
    
    pub fn evolve(&mut self, nb_step : usize)
    { // function of evolution of the GRN, with change of all the proteins concentration, with respect of the nb steps taken
        for _ in 0..nb_step 
        {
            let mut next_proteins: Vec<GrnProtein> = Vec::new();
            // Compute the new proteins
            for j in 0..self.proteins.len() {
                if self.proteins[j].type_ == grn_protein::INPUT_PROTEIN {
                    next_proteins.push(self.proteins[j].clone());
                } else {
                    let mut enhance = 0.0;
                    let mut inhibit = 0.0;

                    for k in 0..self.proteins.len() {
                        if self.proteins[k].type_ != grn_protein::OUTPUT_PROTEIN {
                            enhance += self.proteins[k].concentration * self.enhance_matching[j][k];
                            inhibit += self.proteins[k].concentration * self.inhibit_matching[j][k];
                        }
                    }

                    let new_concentration = (self.proteins[j].concentration + self.delta / self.proteins.len() as f64 * (enhance - inhibit)).max(0.0);
                    next_proteins.push(GrnProtein {
                        id: self.proteins[j].id,
                        type_: self.proteins[j].type_,
                        concentration: new_concentration,
                        enhancer: self.proteins[j].enhancer,
                        inhibiter: self.proteins[j].inhibiter,
                    });
                }
            }

            // total concentration computation
            let mut sum_concentration = 0.0;
            for protein in &next_proteins {
                if protein.type_ != grn_protein::INPUT_PROTEIN {
                    sum_concentration += protein.concentration;
                }
            }

            // concentration normalization 
            if sum_concentration != 0.0 {
                for protein in &mut next_proteins {
                    if protein.type_ != grn_protein::INPUT_PROTEIN {
                        protein.concentration /= sum_concentration;
                    }
                }
            }

            self.proteins = next_proteins;
            self.current_step += 1;
        }
    }

    pub fn distance_to(&self, g : &GrnModel, compare_dynamics_coeff : bool, prot_coef : f64, enh_coef : f64 , 
        inh_coef : f64, beta_max : f64, beta_min : f64, delta_max : f64, delta_min : f64 ) -> f64
    {
        // if g.is_none() {
        //     eprintln!("(GRNModel.distanceTo) Genome is null!");
        // }
        let mut distance = 0.0;
        let gs : &GrnModel;
        let gl : &GrnModel;
        
        
        if self.get_size() < g.get_size() {
            gs = self;
            gl = g;
        } else {
            gs = g;
            gl = self;
        }
        
        // comparing all types
        for gi1 in &gl.proteins {
            let mut min_d = f64::MAX;
            for gi2 in &gs.proteins {
                let d = gi1.distance_to(gi2, prot_coef, enh_coef, inh_coef);
                if gi1.type_ == gi2.type_ && d < min_d {
                    min_d = d;
                }
            }
            distance += min_d;
        }
        
        if compare_dynamics_coeff {
            // take beta and delta to the distance calculation
            distance += (self.beta - g.beta).abs() / (beta_max - beta_min);
            distance += (self.delta - g.delta).abs() / (delta_max - delta_min);
            distance / (gl.get_size() as f64 + 2.0)
        } else {
            distance / gl.get_size() as f64
        }
    }

    // getter and setter 
    pub fn  get_proteins(&self) -> &Vec<GrnProtein> {
        &self.proteins
    }

    pub fn reset(&mut self)
    {
        let prot_size = self.proteins.len();
        for p in &mut self.proteins
        { p.concentration = 1./(prot_size as f64);}
    }    

    fn get_size(&self) -> usize
    {
        self.proteins.len()
    }



}
        