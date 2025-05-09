use core::str;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use rand::prelude::*;

pub trait Graph {
    // fn new() -> Self;
    // fn new_with_size(vertex_ct: usize) -> Self;

    fn resize(&mut self, vertex_ct: usize);
    fn clear(&mut self);

    fn vertex_ct(&self) -> usize;
    fn edge_ct(&self) -> usize;

    fn adjacent(&self, u: usize, v: usize) -> bool;
    fn edges(&self) -> impl Iterator<Item = (usize, usize)>;
    fn vertices(&self) -> impl Iterator<Item = usize>;
    fn neighbors(&self, v: usize) -> impl Iterator<Item = usize>;
    fn degree(&self, v: usize) -> usize {
        self.neighbors(v).count()
    }

    fn add_vertex(&mut self) -> usize;
    fn add_edge(&mut self, u: usize, v: usize);
    fn merge(&mut self, g2: &Self);

    fn remove_edge(&mut self, u: usize, v: usize);

    fn read_dimacs(&mut self, filename: &str) -> Result<(), std::io::Error> {
        let input_file = File::open(filename)?;
        let mut input_reader = BufReader::new(input_file);
        let mut line = String::new();
        let mut line_ctr = 1;

        input_reader.read_line(&mut line)?;
        let mut fields = line.split_whitespace().map(|s| s.parse::<i64>());
        let vertex_ct = fields.next().expect("Could not parse vertex count.").unwrap() as usize;
        let _edge_ct = fields.next().expect("Could not parse edge count.").unwrap() as usize;

        self.clear();
        self.resize(vertex_ct);

        line.clear();
        while input_reader.read_line(&mut line)? > 0 {
            line_ctr += 1;
            let mut fields = line.split_whitespace().map(|s| s.parse::<i64>().expect(format!("Could not parse value on line {}", line_ctr).as_str()));
            
            let u = fields.next().expect(format!("Not enough fields on line {}", line_ctr).as_str()) as usize;
            let v = fields.next().expect(format!("Not enough fields on line {}", line_ctr).as_str()) as usize;

            self.add_edge(u, v);

            line.clear();
        }

        Ok(())
    }

    fn write_dimacs(&self, filename: &str) -> Result<(), std::io::Error> {
        let output_file = File::create(filename)?;
        let mut output_writer = BufWriter::new(output_file);

        writeln!(&mut output_writer, "{}\t{}", self.vertex_ct(), self.edge_ct())?;

        for (u, v) in self.edges() {
            writeln!(&mut output_writer, "{}\t{}", u, v)?;
        }

        Ok(())
    }

    fn random(&mut self, vertex_ct: usize, density: f64) {
        assert!(0. <= density && density <= 1.);
        let mut rng = rand::rng();
        let num_edges = ((vertex_ct * (vertex_ct - 1) / 2) as f64 * density).ceil() as usize;
        let mut edges_added = 0usize;
        self.clear();
        self.resize(vertex_ct);

        while edges_added < num_edges {
            let u = rng.next_u64() as usize % vertex_ct;
            let v = rng.next_u64() as usize % vertex_ct;

            if u != v && !self.adjacent(u, v) {
                self.add_edge(u, v);
                edges_added += 1;
            }
        }
    }

    fn random_regularish(&mut self, vertex_ct: usize, degree: usize) {
        let mut rng = rand::rng();
        self.clear();
        self.resize(vertex_ct);

        for v in 0..vertex_ct {
            while self.degree(v) < degree {
                let mut u = rng.next_u64() as usize % vertex_ct;
                while self.degree(u) > degree || u == v {
                    u = rng.next_u64() as usize % vertex_ct;
                }
                self.add_edge(u, v);
            }
        }
    }


}

pub struct AdjacencyMatrix {
    am: Vec<u8>,
    vertex_ct: usize,
}

impl AdjacencyMatrix {
    pub fn new() -> AdjacencyMatrix {
        let am = Vec::new();
        let vertex_ct = 0;

        AdjacencyMatrix {am, vertex_ct}
    }

    fn set_edge(&mut self, u: usize, v: usize, value: u8) {
        // I have to do this because Rust won't let me subtract with overflow even when I am multiplying the result by 0
        let u = u as i64;
        let v = v as i64;
        if u > v {
            self.am[(u * (u - 1) / 2 + v) as usize] = value;
        }
        else if v > u {
            self.am[(v * (v - 1) / 2 + u) as usize] = value;
        }
        else {
            panic!("Can't add self loops in simple graph.");
        }
    }

    pub fn print_matrix(&self) {
        for u in 0..self.vertex_ct {
            for v in 0..u {
                print!("{}", if self.adjacent(u, v) { 1 } else { 0 });
            }
            println!()
        }
    }
}

impl Graph for AdjacencyMatrix {
    fn resize(&mut self, vertex_ct: usize) {
        self.vertex_ct = vertex_ct;
        self.am.resize(vertex_ct * (vertex_ct - 1) / 2, 0);
    }
    
    fn clear(&mut self) {
        self.am.clear();
    }

    fn vertex_ct(&self) -> usize {
        self.vertex_ct
    }

    fn edge_ct(&self) -> usize{
        match self.am.iter().map(|reference| *reference as usize).reduce(|total, element| (total + element)) {
            None => 0,
            Some(edge_ct) => edge_ct
        }
    }

    fn adjacent(&self, u: usize, v: usize) -> bool {
        let u = u as i64;
        let v = v as i64;
        if u > v {
            self.am[(u * (u - 1) / 2 + v) as usize] == 1
        }
        else if v > u {
            self.am[(v * (v - 1) / 2 + u) as usize] == 1
        }
        else {
            false
        }
    }

    fn edges(&self) -> impl Iterator<Item = (usize, usize)> {
        (0..self.vertex_ct).map(|u| (0..u).map(move |v| (u, v))).flatten().filter(|(u, v)| self.adjacent(*u, *v))
    }

    fn vertices(&self) -> impl Iterator<Item = usize> {
        0..self.vertex_ct
    }

    fn neighbors(&self, v: usize) -> impl Iterator<Item = usize> {
        (0..self.vertex_ct).filter(move |u| self.adjacent(*u, v))
    }

    // fn degree(&self, v: usize) -> usize {
    //     let mut degree = 0usize;

    //     for u in 0..v {
    //         degree += self.am[v * (v - 1) / 2 + u] as usize;
    //     } 

    //     for u in (v + 1)..self.vertex_ct {
    //         degree += self.am[u * (u - 1) / 2 + v] as usize;
    //     }

    //     degree
    // }

    fn add_vertex(&mut self) -> usize {
        self.vertex_ct += 1;
        self.am.resize(self.vertex_ct * (self.vertex_ct - 1), 0);
        return self.vertex_ct - 1;
    }

    fn add_edge(&mut self, u: usize, v: usize) {
        self.set_edge(u, v, 1);
    }

    fn remove_edge(&mut self, u: usize, v: usize) {
        self.set_edge(u, v, 0);
    }

    fn merge(&mut self, g: &Self) {
        let am_len = self.am.len();
        let vertex_ct = self.vertex_ct;
        self.resize(self.vertex_ct + g.vertex_ct);

        // Copy over matrix
        // g.am[..self.am.len()].copy_from_slice(&self.am[..]);

        for v in 1..g.vertex_ct {
            self.am[(am_len + (v + 1) * vertex_ct + v * (v - 1) / 2)..(am_len + (v + 1) * vertex_ct + v * (v + 1) / 2 )]
                .copy_from_slice(&g.am[(v * (v - 1) / 2)..((v + 1) * v / 2)]);
        }

    }
}
