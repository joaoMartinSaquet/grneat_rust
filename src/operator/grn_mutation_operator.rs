use crate::evolver::grn_genome::GrnGenome;


trait grn_mutation_operator {
    fn get_probability(&self) -> f64;
    fn set_probability(&self);
    fn mutate_by_modifying(&self, a_genome : GrnGenome) -> bool;
    fn clone_and_mutate(&self, a_genome : GrnGenome) -> GrnGenome{
        let offspring :GrnGenome  = a_genome.clone();
        if self.mutate_by_modifying(offspring.clone())
        { 
            offspring
        }else {
            a_genome
        }
        
    }
}

