use std::error::Error;
use cfx::Big;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();

    // Expect: modexp <base> <exp> <mod>
    if args.len() != 4 {
        eprintln!("usage: modexp <base> <exp> <mod>");
        eprintln!("example: modexp 4 13 497");
        std::process::exit(2);
    }

    let base = Big::from_dec(&args[1]).ok_or("invalid base")?;
    let exp  = Big::from_dec(&args[2]).ok_or("invalid exp")?;
    let modu = Big::from_dec(&args[3]).ok_or("invalid modulus")?;

    let out = base.mod_exp(&exp, &modu);

    println!("{}^{} mod {} = {}", base, exp, modu, out);
    Ok(())
}
