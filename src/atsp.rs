use std::{collections::{HashMap, HashSet}, mem::swap, ops::Index};

use nalgebra::Vector2;

use crate::graph::{AdjacencyMatrix, Graph, GraphAL};

const C0: f32 = 0.2;

fn compute_r0(points: &Vec<Vector2<f32>>) -> f32 {
    5. * points.iter()
        .map(|point| point.norm())
        .max_by(|a, b| a.partial_cmp(b).expect("Bad values in point list"))
        .expect("Point list empty")
}

fn next_n_k(points: &Vec<Vector2<f32>>, original_net: &Vec<usize>, remaining_points: &HashSet<usize>, r0: f32) -> i32 {
    // d: The maximum distance from a point not in the net to the closest one in the net
    let mut d = 0.0f32;
    
    for new_point in remaining_points {
        let mut min_dist = f32::MAX;
        for og_point in original_net {
            let dist = (points[*og_point] - points[*new_point]).norm();
            if dist < min_dist {
                min_dist = dist;
            }
        }

        if min_dist > d {
            d = min_dist;
        }
    }

    // Choose the smallest k such that (2 ^ -k) * epsilon < d
    (-(d / r0).log2()).ceil() as i32
}

fn next_net(points: &Vec<Vector2<f32>>, original_net: &Vec<usize>, remaining_points: &mut HashSet<usize>, r0: f32, n: i32) -> Vec<usize> {
    let mut new_net = original_net.clone();
    let r = r0 / ((1 << n) as f32);
    let r2 = r * r;
    let mut outside_points = remaining_points.iter().filter(|rem_point| {
        for og_point in original_net {
            if (points[*og_point] - points[**rem_point]).norm_squared() < r * r {
                return false;
            }
        }
        true
    }).map(|i| *i).collect::<Vec<usize>>();

    while outside_points.len() > 0 {
        let new_net_point = outside_points.pop().unwrap();
        new_net.push(new_net_point);
        remaining_points.remove(&new_net_point);
        
        outside_points = outside_points.iter().filter(|i| {
            (points[**i] - points[new_net_point]).norm_squared() >= r2
        }).map(|v| *v).collect::<Vec<usize>>();
    }

    new_net
}

#[derive(Default)]
struct Cylinder {
    x0: Vector2<f32>,
    x1: Vector2<f32>,
    width: f32,
}

fn max_distance_from_line(points: &Vec<Vector2<f32>>, x0: Vector2<f32>, x1: Vector2<f32>, in_ball: &HashSet<usize>) -> f32 {

    match in_ball.iter().map(|index| ((points[*index] - x0) - ((points[*index] - x0).dot(&(x1 - x0)) / (x1 - x0).norm_squared()) * (x1 - x0)).norm()).max_by(|x, y| x.partial_cmp(y).unwrap()) {
        Some(max) => max,
        None => 0.0
    }
  
}

fn get_points_in_ball(points: &Vec<Vector2<f32>>, net: &Vec<usize>, center: usize, radius: f32) -> HashSet<usize> {
    net.iter().map(|index_ref| *index_ref).filter(|index| (points[*index] - points[center]).norm_squared() < radius * radius).collect()
}

fn find_thinnest_cylinder(points: &Vec<Vector2<f32>>, net: &Vec<usize>, center: usize, radius: f32) -> Cylinder {
    let side_offsets = vec![radius * Vector2::<f32>::x(), radius * Vector2::<f32>::y(), -radius * Vector2::<f32>::x(), -radius * Vector2::<f32>::y()];
    let l = (40. * C0).ceil() as i32;
    let mut cylinder_width = f32::MAX;
    let points_in_ball = get_points_in_ball(points, net, center, radius);
    let mut cylinder = Cylinder::default();

    // Atrocious nesting
    for side1 in 0..side_offsets.len() {
        for side2 in 0..side1 {
            for i in (-2 * l + 1)..(2 * l) {
                let x0 = side_offsets[side1] + Vector2::<f32>::new(side_offsets[side1].y, -side_offsets[side1].x) * radius * ((i as f32) / ((2 * l) as f32));
                for j in (-2 * l + 1)..(2 * l) {
                    let x1 = side_offsets[side2] + Vector2::<f32>::new(side_offsets[side2].y, -side_offsets[side2].x) * radius * ((j as f32) / ((2 * l) as f32));

                    let distance = max_distance_from_line(points, x0, x1, &points_in_ball);

                    if distance > cylinder_width {
                        cylinder_width = distance;
                        cylinder.width = cylinder_width;
                        cylinder.x0 = x0;
                        cylinder.x1 = x1;
                    }
                }
            }
        }
    }

    cylinder
}


fn get_bounding_cylinders(points: &Vec<Vector2<f32>>, net :&Vec<usize>, next_net: &Vec<usize>, radius: f32) -> HashMap<usize, Cylinder> {
    let mut cylinder_map = HashMap::new();

    for index in net {
        cylinder_map.insert(*index, find_thinnest_cylinder(points, next_net, *index, radius));
    }

    cylinder_map
}

fn flatness(cylinder: &Cylinder, radius: f32) -> f32{
    cylinder.width / radius
}

fn component(v: &Vector2<f32>, along: &Vector2<f32>) -> f32 {
    v.dot(along) / along.norm()
}

// Step k + 1
// Take a look at what fields are available here and try to implement each part so that you are using somewhat analogous structure
// The paper calls for remaining points from each net individually but I don't think they are used so I only keep track of the 
// Remaining points after the last net
fn astp_step_k_plus_1(
    points: &Vec<Vector2<f32>>,
    graph: &AdjacencyMatrix,
    net_k: &Vec<usize>,
    net_k_plus_1: &Vec<usize>,
    flatness_map_k: &HashMap<usize, Cylinder>,
    n_k: i32,
    n_k_plus_1: i32,
    remaining_points: &mut HashSet<usize>,
    r0: f32,
) {
    let n_k_plus_2 = next_n_k(points, net_k_plus_1, remaining_points, r0);
    let net_k_plus_2 = next_net(points, net_k_plus_1, remaining_points, r0, n_k_plus_2);
    let mut next_graph = AdjacencyMatrix::new();
    next_graph.resize(net_k_plus_1.len());

    if remaining_points.is_empty() {
        return;
    }

    // 5.2.1 Families N_{k+1} and F_{k+1}
    let flatness_map_k_plus_1 = get_bounding_cylinders(points, net_k_plus_1, &net_k_plus_2, r0 * 2.0f32.powi(-(n_k_plus_1 as i32)));

    // 5.2.2 Edges coming from E_k
    let max_edge_length  = C0 * 2.0f32.powi(-n_k_plus_1 - 1) * r0;
    for (u, v) in graph.edges() {
        let flatness_u = flatness_map_k.get(&u).unwrap().width / (r0 * 2.0f32.powi(-n_k_plus_1));
        let flatness_v = flatness_map_k.get(&v).unwrap().width / (r0 * 2.0f32.powi(-n_k_plus_1));

        if (points[u] - points[v]).norm() >= max_edge_length || (flatness_u >= 1. / 16. && flatness_v >= 1. / 16.) {
            next_graph.add_edge(u, v);
        }
        else {
            let mut u = u;
            let mut v = v;
            if flatness_u < 1. / 16. {
                swap(&mut u, &mut v);
            }

            let mut points_along_edge: Vec<(usize, f32)> = get_points_in_ball(points, net_k, u, 2.0f32.powi(-n_k)).union(&get_points_in_ball(points, net_k, v, 2.0f32.powi(-n_k))).into_iter()
                .map(|index| (*index, component(&(points[*index] - points[u]), &(points[v] - points[u])))).collect();
            points_along_edge.sort_by(|(_i1, comp1), (_i2, comp2)| comp1.partial_cmp(comp2).unwrap());

            let index_of_first_past_u = match points_along_edge.binary_search_by(|(_i, comp1)| comp1.partial_cmp(&0.0f32).unwrap()) { Ok(index) => index + 1, Err(index) => index};
            for i in index_of_first_past_u..points_along_edge.len() {
                let (index, _component) = points_along_edge[i];
                let (prev_index, _prev_component) = points_along_edge[i - 1];

                next_graph.add_edge(index, prev_index);

                if index == v {
                    break;
                }
            }
        }
    }
}


// This needs to be rewritten
pub fn atsp(points: &Vec<Vector2<f32>>) -> Vec<i32> {
    let r0 = compute_r0(&points);
    let mut n_k = 1; // n_1 from the paper
    let mut net = vec![0]; // V_1 from the paper
    let mut remaining_points: HashSet<usize> = (1..points.len()).collect();
    let mut i = 1;
    let mut g = AdjacencyMatrix::new();
    g.resize(points.len());

    println!("X_{}: {:?}", i, net);
    println!("n_{}: {}", i, n_k);
    println!("{} points remain", remaining_points.len());

    while remaining_points.len() > 0 {

        let n_k_plus_1 = next_n_k(points, &net, &remaining_points, r0);
        let next_net = next_net(points, &net, &mut remaining_points, r0, n_k_plus_1);
        let bounding_cylinders = get_bounding_cylinders(points, &net, &next_net, r0 * 2.0f32.powi(-(n_k_plus_1 as i32)));
        

        net = next_net;
        n_k = n_k_plus_1;
        i += 1;

        println!("X_{}: {:?}", i, net);
        println!("n_{}: {}", i, n_k);
        println!("{} points remain", remaining_points.len());
    }

    // find_thinnest_cylinder(Vector2::<f32>::zeros(), points, 1.0);

    Vec::new()
}

