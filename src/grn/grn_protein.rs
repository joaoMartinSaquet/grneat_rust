// Proteins constant
pub const INPUT_PROTEIN : i32 = 0;
pub const REGULATORY_PROTEIN: i32 = 2;
pub const OUTPUT_PROTEIN: i32  = 1;
pub const IDSIZE : i32 = 32;


// string indexes 
const ID_IND : usize = 0;
const ENH_IND : usize = 1;
const INH_IND : usize = 2;
const CONCENTRATION_IND : usize = 3;
const TYPE_IND : usize = 4;


#[derive(Debug, Clone, Copy)]
pub struct GrnProtein
{
    pub id: i32,
    pub concentration : f64,
    pub enhancer : i32,
    pub inhibiter : i32,
    pub type_ : i32,

}


impl GrnProtein {
    pub fn new() -> GrnProtein
    {
        GrnProtein {
            id : 0, 
            concentration : 0.0,
            enhancer : 0,
            inhibiter : 0,
            type_ : 0,
        }
    }

    pub fn from_string(prot : String) -> GrnProtein
    {
        let cleaned_prot = prot.replace(&['[', ']'][..], "");
        let splitted_values : Vec<&str> = cleaned_prot.split(',').collect();
        
        let id_ = splitted_values[ID_IND].parse::<i32>().unwrap();
        let enhancer_ = splitted_values[ENH_IND].parse::<i32>().unwrap();
        let inhibiter_ = splitted_values[INH_IND].parse::<i32>().unwrap();
        let concentration_ = splitted_values[CONCENTRATION_IND].parse::<f64>().unwrap();
        let type__ = splitted_values[TYPE_IND].parse::<i32>().unwrap();
        
        GrnProtein {
            id : id_,
            concentration : concentration_,
            enhancer : enhancer_,
            inhibiter : inhibiter_,
            type_ : type__,
        }
    } 

    pub fn distance_to(&self, g: &GrnProtein, prot_coef: f64, enh_coef: f64, inh_coef: f64) -> f64 {
        // if g.is_none() {
        //     eprintln!("(GRNProtein.distanceTo) g is null!!!");
        // }
        (prot_coef * (self.id - g.id).abs() as f64 +
         enh_coef * (self.enhancer - g.enhancer).abs() as f64 +
         inh_coef * (self.inhibiter - g.inhibiter).abs() as f64) / IDSIZE as f64
    }
}