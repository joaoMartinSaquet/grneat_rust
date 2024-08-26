use crate::{evolver::grn_gene::GrnGene, grn::grn_protein};
use std::collections::{hash_map, HashMap};
use rand::*;

pub struct GrnGenome {
    input_genes : HashMap<i32, GrnGene>,
    output_genes : HashMap<i32, GrnGene>,
    regulatory_genes : HashMap<i64, GrnGene>,
    all_genes : HashMap<i64, GrnGene>,

    beta : f64,
	delta : f64,
	beta_max : f64,
	beta_min : f64,
	delta_max : f64,
	delta_min : f64,

    pub parent_best_fit: f64,
    adjusted_fit : f64,
    last_fit : f64,
    has_been_evaluated : bool,
}

impl GrnGenome {

    // constructors 
    pub fn new() -> GrnGenome {
        GrnGenome {
            input_genes : HashMap::<i32, GrnGene>::new(),
            output_genes : HashMap::<i32, GrnGene>::new(),
            regulatory_genes : HashMap::<i64, GrnGene>::new(),
            all_genes : HashMap::<i64, GrnGene>::new(),
            beta : 1.,
            delta : 1.,
            beta_max : 2.0,
            beta_min : 0.5,
            delta_max : 2.,
            delta_min : 0.5,
            parent_best_fit : 0.,
            adjusted_fit: 0.,
            last_fit : 0.,
            has_been_evaluated : false,
        }
    }

    pub fn with_vec(genes : Vec<GrnGene>, beta : f64, delta : f64) -> GrnGenome {
        let mut genome = Self::new();
        let mut iter_gene = genes.iter();
        while iter_gene.len()>0 {
            if let Some(g) = iter_gene.next()
            {   
                genome.add_gene(g.clone());
            }
        }
        genome.beta = beta;
        genome.delta = delta;
        genome
    }

    // methods
    pub fn add_gene(&mut self, a_gene : GrnGene){
        let g = a_gene.clone();
        match g.get_prot_type(){
            grn_protein::INPUT_PROTEIN      => self.input_genes.insert(g.get_connect_to(), g.clone()),
            grn_protein::REGULATORY_PROTEIN => self.input_genes.insert(g.get_id(), g.clone()),
            grn_protein::OUTPUT_PROTEIN     => self.input_genes.insert(g.get_connect_to(), g.clone()),
            _                               => unimplemented!("should not exist unknown protein type"),
        };
        self.all_genes.insert(g.get_id() as i64, g.clone());
        self.has_been_evaluated = true;
    }

    pub fn compare_to(&self, other: &GrnGenome) -> i32 {
        let mut distance = 0;

        // Comparing input genes
        for gi1 in self.input_genes.values() {
            for gi2 in other.input_genes.values() {
                if gi1.compare_to(gi2) != 0 {
                    distance += 1;
                }
            }
        }

        // Comparing output genes
        for go1 in self.output_genes.values() {
            for go2 in other.output_genes.values() {
                if go1.compare_to(go2) != 0 {
                    distance += 1;
                }
            }
        }

        // Comparing regulatory genes
        for gr1 in self.regulatory_genes.values() {
            for gr2 in other.regulatory_genes.values() {
                if gr1.compare_to(gr2) != 0 {
                    distance += 1;
                }
            }
        }

        distance
    }

    pub fn distance_to(&self, g: &GrnGenome, compare_dynamics_coeff: bool) -> f64 {
        // if g.is_null() {
        //     eprintln!("(GRNGenome.distanceTo) Genome is null!");
        // }
        
        let mut distance = 0.0;
        let (gs, gl) = if self.size() < g.size() {
            (self, g)
        } else {
            (g, self)
        };

        // Comparing inputs
        for gi1 in gl.input_genes.values() {
            if let Some(gi2) = gs.input_genes.get(&gi1.get_connect_to()) {
                distance += gi1.distance_to(gi2);
            } else {
                eprintln!("(GRNGenome.distanceTo) gi2 is null! gi1.connectTo={}", gi1.get_connect_to());
                eprintln!("{}", gs.to_string());
            }
        }

        // Comparing outputs
        for go1 in gl.output_genes.values() {
            if let Some(go2) = gs.output_genes.get(&go1.get_connect_to()) {
                distance += go1.distance_to(go2);
            }
        }

        // Comparing regulatory genes
        for gr1 in gl.regulatory_genes.values() {
            let mut min_dist = f64::MAX;
            for gr2 in gs.regulatory_genes.values() {
                min_dist = min_dist.min(gr1.distance_to(gr2));
            }
            distance += min_dist;
        }

        if compare_dynamics_coeff {
            // Take beta and delta into the distance calculation
            distance += (self.beta - g.beta).abs() / (self.beta_max - self.beta_min);
            distance += (self.delta - g.delta).abs() / (self.delta_max - self.delta_min);
            distance / (gl.size() as f64 + 2.0)
        } else {
            distance / gl.size() as f64
        }
    }

    pub fn size(&self) -> usize {
        self.all_genes.len()
    }

    pub fn to_string(&self) -> String {
        let mut res = String::from("{");

        // Concatenate all gene strings
        for g in self.all_genes.values() {
            res.push_str(&g.to_string());
        }

        // Append beta and delta
        res.push_str(&format!(" ; {:.2} ; {:.2}", self.beta, self.delta));

        res
    }

    pub fn contains_gene(&self, genes : GrnGene) -> GrnGene {
        let gene_id = genes.get_id() as i64;
        self.all_genes.get(&gene_id).unwrap().clone()
    }

    pub fn contains_gene_id(&self, genes_id : i64) -> GrnGene {
        self.all_genes.get(&genes_id).unwrap().clone()
    }


    // getter and setter 
    pub fn get_input_genes(&self) -> &HashMap<i32, GrnGene> {
        &self.input_genes
    }

    pub fn get_output_genes(&self) -> &HashMap<i32, GrnGene> {
        &self.output_genes
    }

    pub fn get_regulatory_genes(&self) -> &HashMap<i64, GrnGene> {
        &self.regulatory_genes
    }

    pub fn get_all_genes(&self) -> &HashMap<i64, GrnGene> {
        &self.all_genes
    }

    pub fn get_beta(&self) -> f64 {
        self.beta
    }

    pub fn get_delta(&self) -> f64 {
        self.delta
    }

    pub fn get_beta_max(&self) -> f64 {
        self.beta_max
    }

    pub fn get_beta_min(&self) -> f64 {
        self.beta_min
    }

    pub fn get_delta_max(&self) -> f64 {
        self.delta_max
    }

    pub fn get_delta_min(&self) -> f64 {
        self.delta_min
    }

    pub fn get_parent_best_fit(&self) -> f64 {
        self.parent_best_fit
    }

    pub fn get_adjusted_fit(&self) -> f64 {
        self.adjusted_fit
    }

    pub fn get_last_fit(&self) -> f64 {
        self.last_fit
    }

    pub fn has_been_evaluated(&self) -> bool {
        self.has_been_evaluated
    }

    pub fn set_input_genes(&mut self, input_genes: HashMap<i32, GrnGene>) {
        self.input_genes = input_genes;
    }

    pub fn set_output_genes(&mut self, output_genes: HashMap<i32, GrnGene>) {
        self.output_genes = output_genes;
    }

    pub fn set_regulatory_genes(&mut self, regulatory_genes: HashMap<i64, GrnGene>) {
        self.regulatory_genes = regulatory_genes;
    }

    pub fn set_all_genes(&mut self, all_genes: HashMap<i64, GrnGene>) {
        self.all_genes = all_genes;
    }

    pub fn set_beta(&mut self, beta: f64) {
        self.beta = beta;
    }

    pub fn set_delta(&mut self, delta: f64) {
        self.delta = delta;
    }

    pub fn set_beta_max(&mut self, beta_max: f64) {
        self.beta_max = beta_max;
    }

    pub fn set_beta_min(&mut self, beta_min: f64) {
        self.beta_min = beta_min;
    }

    pub fn set_delta_max(&mut self, delta_max: f64) {
        self.delta_max = delta_max;
    }

    pub fn set_delta_min(&mut self, delta_min: f64) {
        self.delta_min = delta_min;
    }

    pub fn set_parent_best_fit(&mut self, parent_best_fit: f64) {
        self.parent_best_fit = parent_best_fit;
    }

    pub fn set_adjusted_fit(&mut self, adjusted_fit: f64) {
        self.adjusted_fit = adjusted_fit;
    }

    pub fn set_last_fit(&mut self, last_fit: f64) {
        self.last_fit = last_fit;
    }

    pub fn set_has_been_evaluated(&mut self, has_been_evaluated: bool) {
        self.has_been_evaluated = has_been_evaluated;
    }

    pub fn get_input_gene_connect_to(&self, sensor_id : i32) -> GrnGene{
		self.input_genes.get(&sensor_id).unwrap().clone()
	}

    pub fn get_output_gene_connect_to(&self, output_id : i32) -> GrnGene {
	    self.output_genes.get(&output_id).unwrap().clone()
	}    

    pub fn remove_all_genes(&mut self) {
        self.has_been_evaluated = false;
        self.all_genes.clear();
    }

    pub fn remove_gene(&mut self, g: &GrnGene) {
        match g.prot_type {
            grn_protein::INPUT_PROTEIN => {
                self.input_genes.remove(&g.connect_to);
            }
            grn_protein::OUTPUT_PROTEIN => {
                self.input_genes.remove(&g.connect_to);
            }
            grn_protein::REGULATORY_PROTEIN => {
                self.input_genes.remove(&g.get_id());
            }
            _ => {
                unimplemented!("should not happen, unknow protein type");
            } 
        }
        self.all_genes.remove(&(g.get_id() as i64));
        self.has_been_evaluated = false;
    }

    pub fn remove_randomly_regulatory_gene(&mut self, rng: &mut impl Rng) -> bool {
        if self.regulatory_genes.is_empty() {
            return false;
        }

        let delete_index = rng.gen_range(0..self.regulatory_genes.len());
        let g = self.regulatory_genes.values().nth(delete_index).unwrap().clone();

        self.remove_gene(&g);
        self.has_been_evaluated = false;
        true
    }
}

// impl Copy for GrnGenome {
//     fn Copy() -> Self
//     {
        
//     }
// }

impl Clone for GrnGenome
{
    fn clone(&self) -> Self
    {
        GrnGenome {
            input_genes        : self.input_genes.clone(),
            output_genes       : self.output_genes.clone(),
            regulatory_genes   : self.regulatory_genes.clone(),
            all_genes          : self.all_genes.clone(),
            beta               : self.beta,
            delta              : self.delta,
            beta_max           : self.beta_max,
            beta_min           : self.beta_min,
            delta_max          : self.delta_max,
            delta_min          : self.delta_min,
            parent_best_fit    : self.parent_best_fit,
            adjusted_fit       : self.adjusted_fit,
            last_fit           : self.last_fit,
            has_been_evaluated : self.has_been_evaluated,
        }
    }
}
impl Eq for GrnGenome {}

impl PartialEq for GrnGenome {
    fn eq(&self, other: &Self) -> bool {
        self.last_fit == other.last_fit
    }
}

impl Ord for GrnGenome {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.last_fit.partial_cmp(&other.last_fit).unwrap()
    }
}

impl PartialOrd for GrnGenome {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
