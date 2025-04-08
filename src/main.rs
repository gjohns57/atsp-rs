use atsp_rs::atsp::atsp;
use nalgebra::Vector2;
use rand::random;

fn main() {
    let mut points = (0..20).map(|_i| {
        Vector2::new(random::<f32>(), random::<f32>())
    }).collect::<Vec<Vector2<f32>>>();
    let max_magnitude = points.iter().map(|point| point.norm()).min_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    points = points.iter().map(|v| v / max_magnitude).collect();

    println!("{:?}", points);
    atsp(&points);

    println!("Hello, world!");
}
