use crate::grn::{grn_protein, grn_protein::GrnProtein};
use rand::*;

static mut next_gene_id : i32 = 0;

// #[derive(Comparable, Debug)]
pub struct GrnGene {

    prot_coef : f64,
    enh_coef : f64,
    inh_coef : f64,
    current_innov_id : i32,
    gene_id : i32,
    disabled : bool,
    prot_id : i32,
    prot_enh : i32,
    prot_inh : i32,
    pub prot_type : i32,
    pub connect_to :  i32,
}

impl GrnGene {
    
    // constructor
    pub fn default() -> GrnGene {
        GrnGene {
            prot_coef : 0.75,
            enh_coef : 0.125,
            inh_coef : 0.125,
            current_innov_id : 0,
            gene_id : 0,
            disabled : false,
            prot_id : 0,
            prot_enh : 0,
            prot_inh : 0,
            prot_type : grn_protein::REGULATORY_PROTEIN,
            connect_to : 0,
        }
    }

    pub fn with_parameters_prot(n_is_disabled: bool, n_prot_id: i32, n_prot_enh: i32, n_prot_inh: i32, 
        n_prot_type: i32, n_connect_to: i32) -> GrnGene {
            unsafe {
            let id = next_gene_id;
            next_gene_id += 1;
            GrnGene {
                prot_coef : 0.75,
                enh_coef : 0.125,
                inh_coef : 0.125,
                current_innov_id : 0,
                gene_id : (id + 1),
                disabled : n_is_disabled,
                prot_id : n_prot_id,
                prot_enh : n_prot_enh,
                prot_inh : n_prot_inh,
                prot_type : n_prot_type,
                connect_to : n_connect_to,
            }
        }
    }

    pub fn with_prot_values(n_prot_id: i32, n_prot_enh: i32, n_prot_inh: i32) -> GrnGene {

        Self::with_parameters_prot(false, n_prot_id, n_prot_enh, n_prot_inh, grn_protein::REGULATORY_PROTEIN, 0)
    }
    
    pub fn clone(&mut self) -> GrnGene {
        Self::with_parameters_prot(self.disabled, self.prot_id, self.prot_enh,
             self.prot_inh, self.prot_type, self.connect_to)
    }
    
    pub fn compare_to(&self, o: &dyn std::any::Any) -> i32 {
        if let Some(g) = o.downcast_ref::<GrnGene>() {
            if g.prot_type == self.prot_type {
                self.prot_id - g.prot_id
            } else {
                i32::MAX
            }
        } else {
            i32::MAX
        }
    }

    pub fn distance_to(&self, g: &GrnGene) -> f64 {
        // if g.is_none() {
        //     eprintln!("(GRNGene.distanceTo) g is null!!!");
        // }
        if g.prot_type == self.prot_type {
            (self.prot_coef * (self.prot_id - g.prot_id).abs() as f64 +
             self.enh_coef * (self.prot_enh - g.prot_enh).abs() as f64 +
             self.inh_coef * (self.prot_inh - g.prot_inh).abs() as f64) / ( grn_protein::IDSIZE as f64)
        } else {
            1.0
        }
    }

    pub fn generate_random_gene(n_prot_type: i32, connect_to: i32) -> GrnGene {
        // , rng: &mut impl Rng
        let mut rng = rand::thread_rng();
        Self::with_parameters_prot(rng.gen_bool(0.5), 
                                  (rng.gen::<f64>() * grn_protein::IDSIZE as f64) as i32, 
                                  (rng.gen::<f64>() * grn_protein::IDSIZE as f64) as i32, 
                                  (rng.gen::<f64>() * grn_protein::IDSIZE as f64) as i32,
                                  n_prot_type, 
                                  connect_to)
    }

    pub fn generate_random_regulatory_gene() -> GrnGene {
        GrnGene::generate_random_gene(grn_protein::REGULATORY_PROTEIN, 0)
    }
    
    pub fn to_string(&self) -> String {
        format!(
            "[{},{},{},{},{},{},{}]",
            self.prot_id,
            self.prot_enh,
            self.prot_inh,
            self.prot_type,
            self.disabled,
            self.connect_to,
            self.gene_id
        )
    }
    
    // getter and setter
    pub fn is_disbaled(&self) -> bool {
        self.disabled
    }

    pub fn enable(&mut self) {
        self.disabled = false;
    }

    pub fn disable(&mut self) {
        self.disabled = true;
    }

    pub fn get_id(&self) -> i32 {
        self.gene_id
    }

    pub fn get_proteins(&self) -> GrnProtein {
        GrnProtein { id : self.prot_id,
            concentration : 0.0, 
            enhancer : self.prot_enh,
            inhibiter : self.prot_inh,
            type_ : self.prot_type    
        }
    }

    pub fn get_connect_to(&self) -> i32 {
        self.connect_to
    }

    pub fn get_prot_id(&self) -> i32 {
        self.prot_id
    }

    pub fn set_prot_id(&mut self, prot_id: i32) {
        self.prot_id = prot_id;
    }

    pub fn get_prot_enh(&self) -> i32 {
        self.prot_enh
    }

    pub fn set_prot_enh(&mut self, prot_enh: i32) {
        self.prot_enh = prot_enh;
    }

    pub fn get_prot_inh(&self) -> i32 {
        self.prot_inh
    }

    pub fn set_prot_inh(&mut self, prot_inh: i32) {
        self.prot_inh = prot_inh;
    }

    pub fn get_prot_type(&self) -> i32 {
        self.prot_type
    }

    pub fn set_prot_type(&mut self, prot_type: i32) {
        self.prot_type = prot_type;
    }

    pub fn set_disabled(&mut self, is_disabled: bool) {
        self.disabled = is_disabled;
    }
    
}

impl Clone for GrnGene {
    fn clone(&self) -> Self {
            GrnGene {
            prot_coef : self.prot_coef,
            enh_coef : self.enh_coef,
            inh_coef : self.inh_coef,
            current_innov_id : self.current_innov_id,
            gene_id : self.gene_id,
            disabled : self.disabled,
            prot_id : self.prot_id,
            prot_enh : self.prot_enh,
            prot_inh : self.prot_inh,
            prot_type : self.prot_type,
            connect_to : self.connect_to,
        }
    }
}