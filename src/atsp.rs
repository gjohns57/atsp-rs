use std::{
    collections::{HashMap, HashSet, VecDeque},
    mem::swap,
    ops::Index,
};

use nalgebra::{zero, Vector2};

use crate::graph::{AdjacencyMatrix, Graph};

const C0: f32 = 0.2;

fn compute_r0(points: &Vec<Vector2<f32>>) -> f32 {
    5. * points
        .iter()
        .map(|point| point.norm())
        .max_by(|a, b| a.partial_cmp(b).expect("Bad values in point list"))
        .expect("Point list empty")
}

fn next_n_k(
    points: &Vec<Vector2<f32>>,
    original_net: &Vec<usize>,
    remaining_points: &HashSet<usize>,
    r0: f32,
) -> i32 {
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

fn next_net(
    points: &Vec<Vector2<f32>>,
    original_net: &Vec<usize>,
    remaining_points: &mut HashSet<usize>,
    r0: f32,
    n: i32,
) -> Vec<usize> {
    let mut new_net = original_net.clone();
    let r = r0 / ((1 << n) as f32);
    let r2 = r * r;
    let mut outside_points = remaining_points
        .iter()
        .filter(|rem_point| {
            for og_point in original_net {
                if (points[*og_point] - points[**rem_point]).norm_squared() < r * r {
                    return false;
                }
            }
            true
        })
        .map(|i| *i)
        .collect::<Vec<usize>>();

    while outside_points.len() > 0 {
        let new_net_point = outside_points.pop().unwrap();
        new_net.push(new_net_point);
        remaining_points.remove(&new_net_point);

        outside_points = outside_points
            .iter()
            .filter(|i| (points[**i] - points[new_net_point]).norm_squared() >= r2)
            .map(|v| *v)
            .collect::<Vec<usize>>();
    }

    new_net
}

#[derive(Default, Debug)]
struct Cylinder {
    x0: Vector2<f32>,
    x1: Vector2<f32>,
    width: f32,
}

fn max_distance_from_line(
    points: &Vec<Vector2<f32>>,
    x0: Vector2<f32>,
    x1: Vector2<f32>,
    in_ball: &HashSet<usize>,
) -> f32 {
    match in_ball
        .iter()
        .map(|index| {
            ((points[*index] - x0)
                - (points[*index] - x0).dot(&(x1 - x0)) / (x1 - x0).norm_squared() * (x1 - x0))
                .norm()
        })
        .max_by(|x, y| x.partial_cmp(y).unwrap())
    {
        Some(max) => max,
        None => 0.0,
    }
}

fn get_points_in_ball(
    points: &Vec<Vector2<f32>>,
    net: &Vec<usize>,
    center: usize,
    radius: f32,
) -> HashSet<usize> {
    net.iter()
        .map(|index_ref| *index_ref)
        .filter(|index| (points[*index] - points[center]).norm_squared() < radius * radius)
        .collect()
}

fn find_thinnest_cylinder(
    points: &Vec<Vector2<f32>>,
    net: &Vec<usize>,
    center: usize,
    radius: f32,
) -> Cylinder {
    let side_offsets = vec![
        radius * Vector2::<f32>::x(),
        radius * Vector2::<f32>::y(),
        -radius * Vector2::<f32>::x(),
        -radius * Vector2::<f32>::y(),
    ];
    let l = (40. * C0).ceil() as i32;
    let mut cylinder_width = f32::MAX;
    let points_in_ball = get_points_in_ball(points, net, center, radius);
    let mut cylinder = Cylinder::default();
    cylinder.x0 =
        side_offsets[0] + Vector2::<f32>::new(side_offsets[0].y, -side_offsets[0].x) * radius;
    cylinder.x1 =
        side_offsets[2] + Vector2::<f32>::new(side_offsets[2].y, -side_offsets[2].x) * radius;
    cylinder.width = 0.0;

    // Atrocious nesting
    for side1 in 0..side_offsets.len() {
        for side2 in 0..side1 {
            for i in (-2 * l + 1)..(2 * l) {
                let x0 = side_offsets[side1]
                    + Vector2::<f32>::new(side_offsets[side1].y, -side_offsets[side1].x)
                        * radius
                        * ((i as f32) / ((2 * l) as f32));
                for j in (-2 * l + 1)..(2 * l) {
                    let x1 = side_offsets[side2]
                        + Vector2::<f32>::new(side_offsets[side2].y, -side_offsets[side2].x)
                            * radius
                            * ((j as f32) / ((2 * l) as f32));
                    debug_assert_ne!(x0, x1);

                    let distance = max_distance_from_line(points, x0, x1, &points_in_ball);

                    if distance < cylinder_width {
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

fn get_bounding_cylinders(
    points: &Vec<Vector2<f32>>,
    net: &Vec<usize>,
    next_net: &Vec<usize>,
    radius: f32,
) -> HashMap<usize, Cylinder> {
    let mut cylinder_map = HashMap::new();

    for index in net {
        cylinder_map.insert(
            *index,
            find_thinnest_cylinder(points, next_net, *index, radius),
        );
    }

    cylinder_map
}

fn flatness(cylinder: &Cylinder, radius: f32) -> f32 {
    cylinder.width / radius
}

fn component(v: &Vector2<f32>, along: &Vector2<f32>) -> f32 {
    debug_assert_ne!(*along, Vector2::zeros());
    v.dot(along) / along.norm()
}

fn bfs(
    start: usize,
    visited: &mut HashSet<usize>,
    reps: &mut Vec<usize>,
    graph: &mut AdjacencyMatrix,
    used: &mut AdjacencyMatrix,
    set: &HashSet<usize>,
    map: &HashMap<usize, usize>,
    origin: usize,
) {
    let mut q: VecDeque<usize> = VecDeque::new();
    q.push_back(start);

    while !q.is_empty() {
        let node = map[q.front().unwrap()];

        //Pops out the rep if the origin is inside the component
        //So we dont connect origin to its own component
        if node == origin {
            reps.pop();
        }

        //I think flat pairs should be just a line, but
        //kept this just in case
        if visited.contains(&node) {
            continue;
        }

        visited.insert(node);
        q.pop_front();

        for i in set {
            if graph.adjacent(node, map[i]) && !visited.contains(&map[i]) {
                //Checks if the edge is in used
                if used.adjacent(node, map[i]) {
                    //Removes it from the graph if it is
                    graph.remove_edge(node, map[i]);
                } else {
                    //Adds the edge to the used and pushes back the vertex
                    used.add_edge(node, map[i]);
                    q.push_back(map[i]);
                }
            }
        }
    }
}

fn find_reps(
    edges: &mut AdjacencyMatrix,
    used: &mut AdjacencyMatrix,
    set: &HashSet<usize>,
    map: &HashMap<usize, usize>,
    origin: usize,
) -> Vec<usize> {
    let mut reps: Vec<usize> = Vec::new();
    let mut visited: HashSet<usize> = HashSet::new();

    //Goes through each vertex and finds a vertex from each component
    //Will remove edges already used
    for i in set {
        if visited.contains(&map[i]) {
            continue;
        }
        reps.push(map[i]);
        bfs(*i, &mut visited, &mut reps, edges, used, &set, &map, origin);
    }

    reps
}

fn connect(
    origin: usize,
    reps: &Vec<usize>,
    edges: &mut AdjacencyMatrix,
    used: &mut AdjacencyMatrix,
) {
    //Connects the origin to each component rep
    for i in 0..reps.len() {
        edges.add_edge(origin, reps[i]);
        used.add_edge(origin, reps[i]);
    }
}

fn has_edge(used: &AdjacencyMatrix, v: usize) -> bool {
    for i in 0..used.vertex_ct() {
        if used.adjacent(i, v) {
            return true;
        }
    }

    false
}

fn contains_edge(
    graph: &AdjacencyMatrix,
    verts: &HashSet<usize>,
    map: &HashMap<usize, usize>,
) -> HashSet<usize> {
    let mut contains: HashSet<usize> = HashSet::new();

    for v in verts {
        if has_edge(graph, map[v]) {
            contains.insert(*v);
        }
    }

    contains
}

//This function should add edges between flat pairs
fn add_flat_pairs(
    graph: &mut AdjacencyMatrix,
    v: usize,
    points: &Vec<Vector2<f32>>,
    net1: &Vec<usize>,
    net2: &Vec<usize>,
    verts: &mut HashSet<usize>,
    map: &HashMap<usize, usize>,
    epsilon: f32,
    k: i32,
) {
    let xp_in_ball: HashSet<usize> = get_points_in_ball(points, net2, v, epsilon);

    //This feels bad but thinnest cylinder takes a vec
    let xp_in_ball_vec: Vec<usize> = xp_in_ball.into_iter().collect();

    let c: Cylinder = find_thinnest_cylinder(points, &xp_in_ball_vec, v, epsilon);
    let alpha = 2.0f32.powi(k) * c.width / epsilon;

    if alpha < 1.0 / 16. {
        let mut min_cylinder_aligned_component = f32::MAX;

        //Set this to p so that if no flat pair then the first part of the && will fail
        let mut next_point = v;
        for q in xp_in_ball_vec {
            if (points[v] - points[q]).norm() < epsilon
                || (points[q] - points[v]).norm() >= C0 * 2.0f32.powi(-k - 1) * epsilon
            {
                continue;
            }

            let component = component(&(points[q] - points[v]), &(c.x1 - c.x0));
            if component < min_cylinder_aligned_component {
                min_cylinder_aligned_component = component;
                next_point = q;
            }
        }

        if (points[next_point] - points[v]).norm() > epsilon
            && (points[next_point] - points[v]).norm() < 2. * epsilon
        {
            graph.add_edge(map[&v], map[&next_point]);
            verts.insert(v);
            verts.insert(next_point);
        }
    }
}

pub fn non_flat(
    points: &Vec<Vector2<f32>>,
    graph: &mut AdjacencyMatrix, //next_graph
    net_k: &Vec<usize>,
    net_k_plus_1: &Vec<usize>,
    not_flat_k: &Vec<usize>,
    n_k_plus_1: i32,
    n_k: i32,
    r0: f32,
) {
    let mut map: HashMap<usize, usize> = HashMap::new();
    for i in 0..net_k_plus_1.len() {
        map.insert(net_k_plus_1[i], i);
    }

    for u in not_flat_k {
        //Index of all points in net_k_plus_1 that are inside the ball of the current u
        let vp_k_plus_1 =
            get_points_in_ball(points, net_k_plus_1, *u, C0 * 2.0f32.powi(-n_k_plus_1) * r0);

        //Contains edge should return a hashset of the indices that have an edge in next_graph
        //So we take the difference
        let v_k_plus_1: HashSet<usize> = vp_k_plus_1
            .difference(&contains_edge(graph, &vp_k_plus_1, &map))
            .cloned()
            .collect();

        //Goes to next point if its empty
        if v_k_plus_1.is_empty() {
            continue;
        }

        let mut g_k_plus_1 = AdjacencyMatrix::new();
        g_k_plus_1.resize(net_k_plus_1.len());
        let mut in_flat_pair: HashSet<usize> = HashSet::new();

        /* I was not sure what exactly to pass into epsilon and k */
        for v in v_k_plus_1.iter() {
            add_flat_pairs(
                &mut g_k_plus_1,
                *v,
                &points,
                net_k,
                net_k_plus_1,
                &mut in_flat_pair,
                &map,
                r0,
                n_k,
            );
        }

        let vs_k_plus_1: HashSet<usize> = vp_k_plus_1.union(&in_flat_pair).cloned().collect();

        let mut reps = find_reps(&mut g_k_plus_1, graph, &vs_k_plus_1, &map, *u);
        connect(*u, &mut reps, &mut g_k_plus_1, graph);
    }
}

// Step k + 1
// Take a look at what fields are available here and try to implement each part so that you are using somewhat analogous structure
// The paper calls for remaining points from each net individually but I don't think they are used so I only keep track of the
// Remaining points after the last net
fn atsp_step_k_plus_1(
    points: &Vec<Vector2<f32>>,
    graph: &AdjacencyMatrix,
    net_k: &Vec<usize>,
    net_k_plus_1: &Vec<usize>,
    flatness_map_k: &HashMap<usize, Cylinder>,
    n_k: i32,
    n_k_plus_1: i32,
    remaining_points: &mut HashSet<usize>,
    r0: f32,
) -> AdjacencyMatrix {
    println!("{:?}", net_k);
    let mut next_graph = AdjacencyMatrix::new();
    next_graph.resize(net_k_plus_1.len());

    let mut vertex_labels = HashMap::new();
    for i in 0..net_k_plus_1.len() {
        vertex_labels.insert(net_k_plus_1[i], i);
    }

    let r_k = r0 * 2.0f32.powi(-n_k);

    // 5.2.1 Families N_{k+1} and F_{k+1}

    // 5.2.2 Edges coming from E_k
    let max_edge_length = 300. * 2.0f32.powi(-n_k_plus_1 - 1) * r0;
    for (u, v) in graph.edges() {
        println!("Handling edges about edge {{{}, {}}}", net_k[u], net_k[v]);
        let flatness_u = flatness_map_k[&net_k[u]].width / r_k;
        let flatness_v = flatness_map_k[&net_k[v]].width / r_k;

        if (points[net_k[u]] - points[net_k[v]]).norm() >= max_edge_length
            || (flatness_u >= 1. / 16. && flatness_v >= 1. / 16.)
        {
            println!(
                "Keeping edge {} {}",
                vertex_labels[&net_k[u]], vertex_labels[&net_k[v]]
            );
            next_graph.add_edge(vertex_labels[&net_k[u]], vertex_labels[&net_k[v]]);
        } else {
            let mut u = u;
            let mut v = v;
            if flatness_u < 1. / 16. {
                swap(&mut u, &mut v);
            }

            let mut points_along_edge: Vec<(usize, f32)> =
                get_points_in_ball(points, net_k_plus_1, net_k[u], r_k)
                    .union(&get_points_in_ball(points, net_k_plus_1, net_k[v], r_k))
                    .into_iter()
                    .map(|index| {
                        (
                            *index,
                            component(
                                &(points[*index] - points[net_k[u]]),
                                &(points[net_k[v]] - points[net_k[u]]),
                            ),
                        )
                    })
                    .collect();
            points_along_edge
                .sort_by(|(_i1, comp1), (_i2, comp2)| comp1.partial_cmp(comp2).unwrap());

            println!("\tPoints along edge {:?}", points_along_edge);

            // Find the next element after the one with component 0 namely u
            let index_of_first_past_u = match points_along_edge
                .binary_search_by(|(_i, comp1)| comp1.partial_cmp(&0.0).unwrap())
            {
                Ok(index) => index + 1,
                Err(index) => index + 1,
            };

            for i in index_of_first_past_u..points_along_edge.len() {
                let (index, _component) = points_along_edge[i];
                let (prev_index, _prev_component) = points_along_edge[i - 1];

                println!(
                    "\treplacing edge {} {} with {} {}",
                    net_k[u], net_k[v], prev_index, index
                );
                next_graph.add_edge(vertex_labels[&index], vertex_labels[&prev_index]);

                if index == net_k[v] {
                    break;
                }
            }
        }
    }

    // 5.2.3 Edges coming from F_k
    for (u, cyl) in flatness_map_k
        .iter()
        .filter(|(_u, cyl)| cyl.width / r_k < 1.0 / 16.0)
    {
        println!("Adding flat edges about {}", u);
        println!("\tCylinder: {:?}", *cyl);
        let points_in_ball = get_points_in_ball(points, net_k_plus_1, *u, r_k);
        println!("\t{:?}", points_in_ball);

        let net_k_union_ball = points_in_ball.iter().filter(|point| net_k.contains(point));
        let mut do_points_to_right = true;
        let mut do_points_to_left = true;

        for u_prime in net_k_union_ball {
            debug_assert_ne!(cyl.x1, cyl.x0);
            if component(&(points[*u_prime] - points[*u]), &(cyl.x1 - cyl.x0)) < 0. {
                do_points_to_left = false;
            } else if component(&(points[*u_prime] - points[*u]), &(cyl.x1 - cyl.x0)) > 0. {
                do_points_to_right = false;
            }
        }

        if do_points_to_left || do_points_to_right {
            let mut ordered_along_cyl = points_in_ball
                .iter()
                .map(|u1| {
                    (
                        *u1,
                        component(&(points[*u1] - points[*u]), &(cyl.x1 - cyl.x0)),
                    )
                })
                .collect::<Vec<(usize, f32)>>();
            ordered_along_cyl.sort_by(|(_u1, comp1), (_u2, comp2)| {
                comp1
                    .partial_cmp(comp2)
                    .expect("Couldn't compare components")
            });

            let u_index = match ordered_along_cyl.binary_search_by(|(_u1, comp1)| {
                comp1
                    .partial_cmp(&0.0)
                    .expect("Couldn't compare components")
            }) {
                Ok(index) => index,
                Err(index) => index,
            };

            println!("\tordered along cyl {:?}", ordered_along_cyl);

            if do_points_to_left {
                for i in 0..u_index {
                    println!(
                        "\tAdding flat edge {} {}",
                        &ordered_along_cyl[i].0,
                        &ordered_along_cyl[i + 1].0
                    );
                    next_graph.add_edge(
                        vertex_labels[&ordered_along_cyl[i].0],
                        vertex_labels[&ordered_along_cyl[i + 1].0],
                    );
                }
            }
            if do_points_to_right {
                for i in u_index..(ordered_along_cyl.len() - 1) {
                    println!(
                        "\tAdding flat edge {} {}",
                        ordered_along_cyl[i].0,
                        ordered_along_cyl[i + 1].0
                    );

                    next_graph.add_edge(
                        vertex_labels[&ordered_along_cyl[i].0],
                        vertex_labels[&ordered_along_cyl[i + 1].0],
                    );
                }
            }
        }
    }

    // non_flat(points, &mut next_graph, net_k, net_k_plus_1, &flatness_map_k.iter().filter(|(_u, cyl)| cyl.width / r_k < 1.0 / 16.0).map(|tup| *tup.0).collect(), n_k_plus_1, r0);

    if !remaining_points.is_empty() {
        let n_k_plus_2 = next_n_k(points, net_k_plus_1, remaining_points, r0);
        let net_k_plus_2 = next_net(points, net_k_plus_1, remaining_points, r0, n_k_plus_2);
        let flatness_map_k_plus_1 = get_bounding_cylinders(
            points,
            net_k_plus_1,
            &net_k_plus_2,
            r0 * 2.0f32.powi(-(n_k_plus_1 as i32)),
        );

        return atsp_step_k_plus_1(
            points,
            &next_graph,
            net_k_plus_1,
            &net_k_plus_2,
            &flatness_map_k_plus_1,
            n_k_plus_1,
            n_k_plus_2,
            remaining_points,
            r0,
        );
    } else {
        let mut remapped_graph = AdjacencyMatrix::new();
        remapped_graph.resize(next_graph.vertex_ct());

        for (u, v) in next_graph.edges() {
            remapped_graph.add_edge(net_k_plus_1[u], net_k_plus_1[v]);
        }

        return remapped_graph;
    }
}

fn euler_tour(g: &AdjacencyMatrix, traversal: &mut Vec<usize>, visited: &mut Vec<bool>, v: usize) {
    traversal.push(v);
    visited[v] = true;
    for u in g.neighbors(v) {
        if !visited[u] {
            euler_tour(g, traversal, visited, u);
            traversal.push(v);
        }
    }
}

// This needs to be rewritten
pub fn atsp(points: &Vec<Vector2<f32>>) -> Vec<usize> {
    let r0 = compute_r0(&points);
    let n_1 = 1; // n_1 from the paper
    let net = vec![0]; // V_1 from the paper
    let mut remaining_points: HashSet<usize> = (1..points.len()).collect();
    let mut g = AdjacencyMatrix::new();
    g.add_vertex();
    let n_2 = next_n_k(points, &net, &remaining_points, r0);
    let net_2 = next_net(points, &net, &mut remaining_points, r0, n_2);

    let flatness_map_1 = get_bounding_cylinders(points, &net, &net_2, r0);

    let final_graph = atsp_step_k_plus_1(
        points,
        &g,
        &net,
        &net_2,
        &flatness_map_1,
        n_1,
        n_2,
        &mut remaining_points,
        r0,
    );

    let mut traversal = Vec::new();
    let mut visited = Vec::new();
    visited.resize(final_graph.vertex_ct(), false);
    euler_tour(&final_graph, &mut traversal, &mut visited, 0);
    traversal
}