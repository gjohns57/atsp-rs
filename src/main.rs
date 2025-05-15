use atsp_rs::atsp::atsp;
use nalgebra::{vector, Vector2};
use std::io::stdin;
//use rand::random;
fn main() {
    let mut points: Vec<Vector2<f32>> = Vec::new();
    let mut ipt = String::new();
    let mut x: f32;
    let mut y: f32;
    println!("Enter points and then type done!");
    loop {
        stdin().read_line(&mut ipt).expect("unable to read line");

        if ipt.trim() == "done" {
            break;
        }

        let nums: Vec<String> = ipt.trim()
            .split_whitespace()
            .map(String::from)
            .collect();
        if nums.len() != 2 {
            println!("Error, only 2 values per line");
            return;
        }
        x = nums[0].parse::<f32>().expect("unable to read x value");
        y = nums[1].parse::<f32>().expect("unable to read y value");
        points.push(vector![x, y]);
        ipt.clear();
    }
    // let max_magnitude = points.iter().map(|point| point.norm()).max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    // points = points.iter().map(|v| v / max_magnitude).collect();

    println!("{:?}", points);
    let traversal = atsp(&points);

    for u in traversal {
        print!("{} ", u);
    }

    println!();
}