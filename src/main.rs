use atsp_rs::atsp::atsp;
use nalgebra::Vector2;
use rand::random;

fn main() {
    let mut points = (0..5).map(|i| {
        Vector2::new(i as f32, 0.0)
    }).collect::<Vec<Vector2<f32>>>();
    // let max_magnitude = points.iter().map(|point| point.norm()).max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    // points = points.iter().map(|v| v / max_magnitude).collect();

    println!("{:?}", points);
    let traversal = atsp(&points);

    for u in traversal {
        print!("{} ", points[u]);
    }

    println!();
}
