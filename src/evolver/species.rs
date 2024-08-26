use core::f64;
use std::iter::Cloned;
use std::sync::Arc;
use crate::{evaluator, evolver::{grn_gene, grn_genome}};
use grn_genome::GrnGenome;
use crate::evaluator::{grn_genome_evaluator::GrnGenomeEvaluator, grn_sine_evaluator::SineEvaluator};




static mut next_species_id : i32 = 0;


struct Species {
    genomes: Vec<GrnGenome>,
    genomes_are_sorted: bool,
    speciation_threshold: f64,
    representative_genome: Option<GrnGenome>,
    evaluator: Option<SineEvaluator>,
    fitness_sum: f64,
    best_genome: Option<GrnGenome>,
    best_fitness: f64,
    representative_is_first: bool,
    replacement: bool,
    use_dynamics_genome_for_distance: bool,
    species_id: i32, 
}
impl Species {
    fn new() -> Self {
        Species {
            genomes: Vec::<GrnGenome>::new(),
            genomes_are_sorted: true,
            speciation_threshold: 0.3,
            representative_genome: None,
            evaluator:None,
            fitness_sum: 0.0,
            best_genome: None,
            best_fitness: -f64::MAX,
            representative_is_first: true,
            replacement: true,
            use_dynamics_genome_for_distance: true,
            species_id : 0,
        }
    }
    pub fn with_parameters(representative : GrnGenome, speciation_threshold : f64, add_repr : bool, an_evaluator : SineEvaluator) -> Species{
        unsafe {
            let id = next_species_id;
            next_species_id += 1;
        
            let mut s =  Species {
                            genomes: Vec::<GrnGenome>::new(),
                            genomes_are_sorted: true,
                            speciation_threshold: speciation_threshold,
                            representative_genome: Some(representative.clone()),
                            evaluator: Some(an_evaluator),
                            fitness_sum: 0.0,
                            best_genome: None,
                            best_fitness: f64::MIN,
                            representative_is_first: true,
                            replacement: true,
                            use_dynamics_genome_for_distance: true,
                            species_id : id,
                            };
            s.add_genome(representative.clone(), false);
            s
        }
        
    }
	
    pub fn add_genome(&mut self, a_genome : GrnGenome , test_species : bool)  -> bool {
        
        
        if !test_species || self.representative_genome.as_ref().unwrap().distance_to(&a_genome, self.use_dynamics_genome_for_distance)<= self.speciation_threshold {
            if !a_genome.has_been_evaluated() {
                self.evaluator.as_ref().unwrap().evaluate(&a_genome);
            }
            if self.representative_genome.is_none()  {
                self.representative_genome = Some(a_genome.clone());
            }   

            self.genomes.push(a_genome.clone());
            //System.err.println(genomes.size());
            self.genomes_are_sorted=false;
            let best_fit = a_genome.get_last_fit();
            self.fitness_sum += best_fit;
            if (best_fit>self.best_fitness) {
                self.best_genome = Some(a_genome);
                if (!self.representative_is_first) {
                    self.representative_genome=self.best_genome.clone();
                }
            }
                true
            } 
            else {
                false
            }
    } 

	pub fn recompute_fitness(&self) {
	    for k in 0..self.genomes.len() {
            //genomes.set(k, evaluator.evaluate( genomes.get(k) ));
            if let Some(evaluator) = &self.evaluator 
            {   
                evaluator.evaluate( self.genomes.get(k).unwrap() );
            }
	    }
	}   

    pub fn size(&self) -> usize {
        self.genomes.len()
    }

    pub fn remove_worst_genome(&mut self) {
        if !self.genomes_are_sorted {
			self.sort_Genome();
		}
        let last_index = self.genomes.len() - 1;
        self.fitness_sum -= self.genomes.get(last_index).unwrap().get_last_fit();
        self.genomes.remove(last_index);
    }   

    pub fn sort_Genome(&mut self){
        self.genomes.sort();    
        self.genomes_are_sorted = true;
    }

    fn test_sorted(&self) -> bool {
        let mut prev_genome_fit = self.genomes.get(0).unwrap().get_last_fit();
        let mut is_sorted : bool = false;

        for i in 1..self.genomes.len() {
            let current_fit =  self.genomes.get(i).unwrap().get_last_fit();
            is_sorted =  current_fit <= prev_genome_fit; 
            prev_genome_fit = current_fit;
            if !is_sorted { break; }
        }   
        is_sorted
    }

}


impl Clone for Species {
    fn clone(&self) -> Self {
        Species {
            genomes: self.genomes.clone(),
            genomes_are_sorted: self.genomes_are_sorted,
            speciation_threshold: self.speciation_threshold,
            representative_genome: self.representative_genome.clone(),
            evaluator : self.evaluator.clone(),
            // evaluator: self.evaluator.as_ref().map(|evaluator| evaluator),
            fitness_sum: self.fitness_sum,
            best_genome: self.best_genome.clone(),
            best_fitness: self.best_fitness,
            representative_is_first: self.representative_is_first,
            replacement: self.replacement,
            use_dynamics_genome_for_distance: self.use_dynamics_genome_for_distance,
            species_id: self.species_id, 
        }
        
    }
}

// TODO care with static var for next ids, can cause trouble in the implementation