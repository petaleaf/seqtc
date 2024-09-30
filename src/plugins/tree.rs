use crate::core::Plugin;
use clap::{ArgMatches, Command};
use bio::io::fasta;
use std::fs::{File, write};
use std::io::BufReader;
use std::collections::HashMap;

// 定义树节点
#[derive(Debug, Clone)]
struct TreeNode {
    name: String,
    left: Option<Box<TreeNode>>,
    right: Option<Box<TreeNode>>,
    distance: f64,
}

pub struct TreePlugin;

impl Plugin for TreePlugin {
    fn name(&self) -> &'static str {
        "tree"
    }

    fn command(&self) -> Command {
        Command::new("tree")
            .about("构建进化树并输出 Newick 格式")
            .arg(clap::arg!(-m --model <MODEL> "选择建树模型"))
            .arg(
                clap::arg!(-f --fasta <FA_FILE> "输入的FASTQ文件路径")
                .default_value("ML")
            )
            .arg(clap::arg!(-o --out <NEWICK_FILE> "输出的Newick格式文件路径"))
    }

    fn run(&self, matches: &ArgMatches) {
        let model = matches.get_one::<String>("model").expect("请选择需要使用的模型");
        let fa_file = matches.get_one::<String>("fasta").expect("需要输入FASTA文件路径");
        let newick_file = matches.get_one::<String>("out").expect("需要输出Newick文件路径");

        // 读取序列并构建进化树
        match self.build_tree(fa_file ,model, newick_file) {
            Ok(_) => println!("进化树已成功保存到: {}", newick_file),
            Err(err) => eprintln!("构建进化树时发生错误: {}", err),
        }
    }
}

impl TreePlugin {
    fn build_tree(&self, fa_file: &str, model:&str, newick_file: &str) -> Result<(), Box<dyn std::error::Error>> {
        // 读取 FASTQ 序列
        let sequences = read_fasta_sequences(fa_file)?;

        // 生成名称列表
        let names: Vec<String> = sequences.iter().map(|(id, _)| id.clone()).collect();
        let seqs_only: Vec<String> = sequences.iter().map(|(_, seq)| seq.clone()).collect();

        // 计算距离矩阵
        let distance_matrix = compute_distance_matrix(&seqs_only);

        // 构建进化树
        let tree = match model {
            "ML" => build_maximum_likelihood_tree(&seqs_only, &names),
            "NJ" => build_phylogenetic_tree(&distance_matrix, &names),
            _ => return Err(Box::from(format!("不支持的模型类型: {}", model))),
        };
        // 输出为 Newick 格式
        let newick = tree_to_newick(&tree);
        write(newick_file, newick)?;

        Ok(())
    }
}

// 读取 FASTA 文件中的序列
fn read_fasta_sequences(file_path: &str) -> Result<Vec<(String, String)>, Box<dyn std::error::Error>> {
    let reader = BufReader::new(File::open(file_path)?);
    let fasta_reader = fasta::Reader::new(reader);

    let mut sequences: Vec<(String, String)> = Vec::new();
    for result in fasta_reader.records() {
        let record = result?;
        let seq = String::from_utf8_lossy(record.seq()).to_string();
        sequences.push((record.id().to_string(), seq));
    }

    Ok(sequences)
}

// 计算 Hamming 距离
fn hamming_distance(seq1: &str, seq2: &str) -> usize {
    seq1.chars().zip(seq2.chars()).filter(|(a, b)| a != b).count()
}

// 计算距离矩阵
fn compute_distance_matrix(sequences: &[String]) -> Vec<Vec<f64>> {
    let n = sequences.len();
    let mut matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in i..n {
            if i == j {
                matrix[i][j] = 0.0;
            } else {
                let distance = hamming_distance(&sequences[i], &sequences[j]) as f64;
                matrix[i][j] = distance;
                matrix[j][i] = distance; // 距离矩阵是对称的
            }
        }
    }
    matrix
}

// 构建邻接法的进化树
fn build_phylogenetic_tree(distance_matrix: &[Vec<f64>], names: &[String]) -> TreeNode {
    let mut n = distance_matrix.len();
    let mut dist = distance_matrix.to_vec();
    let mut node_names: Vec<String> = names.to_vec();
    let mut tree_nodes: HashMap<String, TreeNode> = HashMap::new();

    // 初始化每个叶节点
    for name in names {
        tree_nodes.insert(
            name.clone(),
            TreeNode {
                name: name.clone(),
                left: None,
                right: None,
                distance: 0.0,
            },
        );
    }

    while n > 2 {
        // 计算邻接法的 Q 矩阵
        let mut q_matrix = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in i + 1..n {
                let row_sum: f64 = dist[i].iter().sum();
                let col_sum: f64 = dist[j].iter().sum();
                q_matrix[i][j] = (n as f64 - 2.0) * dist[i][j] - row_sum - col_sum;
                q_matrix[j][i] = q_matrix[i][j];
            }
        }

        // 寻找最邻近的两个节点
        let mut min_i = 0;
        let mut min_j = 1;
        let mut min_value = q_matrix[0][1];
        for i in 0..n {
            for j in i + 1..n {
                if q_matrix[i][j] < min_value {
                    min_value = q_matrix[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        // 合并两个最近的节点
        let new_node_name = format!("({},{})", node_names[min_i], node_names[min_j]);
        let node_i = tree_nodes.remove(&node_names[min_i]).unwrap();
        let node_j = tree_nodes.remove(&node_names[min_j]).unwrap();

        let distance = dist[min_i][min_j] / 2.0;
        let new_node = TreeNode {
            name: new_node_name.clone(),
            left: Some(Box::new(node_i)),
            right: Some(Box::new(node_j)),
            distance,
        };
        tree_nodes.insert(new_node_name.clone(), new_node);

        // 更新距离矩阵
        let mut new_dist_row = vec![0.0; n - 1];
        for k in 0..n {
            if k != min_i && k != min_j {
                let dist_to_new_node =
                    (dist[min_i][k] + dist[min_j][k] - dist[min_i][min_j]) / 2.0;
                new_dist_row[k.min(n - 2)] = dist_to_new_node;
            }
        }

        // 删除行和列
        dist.remove(min_j);
        dist.remove(min_i);
        for row in &mut dist {
            row.remove(min_j);
            row.remove(min_i);
        }

        dist.push(new_dist_row.clone());
        for i in 0..n - 2 {
            dist[i].push(new_dist_row[i]);
        }

        node_names.remove(min_j);
        node_names.remove(min_i);
        node_names.push(new_node_name);

        n -= 1;
    }

    // 将剩下的两个节点连接为根
    let root_name = format!("({},{})", node_names[0], node_names[1]);
    let node_left = tree_nodes.remove(&node_names[0]).unwrap();
    let node_right = tree_nodes.remove(&node_names[1]).unwrap();
    let root = TreeNode {
        name: root_name.clone(),
        left: Some(Box::new(node_left)),
        right: Some(Box::new(node_right)),
        distance: 0.0,
    };

    root
}


// 实现最大似然法构建进化树
fn build_maximum_likelihood_tree(sequences: &[String], names: &[String]) -> TreeNode {
    let substitution_rate = 0.1; // 假设的替换率，可通过命令行参数调整

    // 将序列存入 HashMap 以便快速访问
    let mut sequence_map: HashMap<String, String> = HashMap::new();
    for (name, seq) in names.iter().zip(sequences.iter()) {
        sequence_map.insert(name.clone(), seq.clone());
    }

    // 构建初始树（贪心法）
    let mut tree = build_initial_tree(&sequence_map);

    // 优化树结构（最大似然法）
    optimize_tree(&mut tree, &sequence_map, substitution_rate);

    tree
}



// Jukes-Cantor 模型替换概率计算
fn jukes_cantor_distance(seq1: &str, seq2: &str) -> f64 {
    let mismatches = seq1.chars().zip(seq2.chars())
        .filter(|(a, b)| a != b)
        .count();
    
    let p = mismatches as f64 / seq1.len() as f64;
    
    // 防止p等于1的情况，保证计算稳定性
    if p >= 0.75 {
        return f64::INFINITY;
    }
    
    // Jukes-Cantor 距离公式
    -0.75 * (1.0 - (4.0 * p / 3.0)).ln()
}


// 定义树结构
impl TreeNode {
    fn new(name: &str) -> Self {
        TreeNode {
            name: name.to_string(),
            left: None,
            right: None,
            distance: 0.0,
        }
    }

    // 设置子节点
    fn set_children(&mut self, left: TreeNode, right: TreeNode, distance: f64) {
        self.left = Some(Box::new(left));
        self.right = Some(Box::new(right));
        self.distance = distance;
    }
}


// 构建初始的树
fn build_initial_tree(sequences: &HashMap<String, String>) -> TreeNode {
    let names: Vec<&String> = sequences.keys().collect();
    let root = TreeNode::new(names[0]);

    let mut current_tree = root;

    for name in names.iter().skip(1) {
        let mut new_node = TreeNode::new(name);
        let distance = jukes_cantor_distance(sequences.get(&current_tree.name).unwrap(), sequences.get(*name).unwrap());
        new_node.set_children(current_tree.clone(), TreeNode::new(name), distance);
        current_tree = new_node;
    }

    current_tree
}

// 贪心法优化树形结构
fn optimize_tree(tree: &mut TreeNode, sequences: &HashMap<String, String>, substitution_rate: f64) {
    let mut best_likelihood = calculate_likelihood(tree, sequences, substitution_rate);

    for _ in 0..100 { // 简单的迭代次数
        // 交换子树，尝试不同的结构
        let swapped_tree = swap_subtrees(tree.clone());
        let likelihood = calculate_likelihood(&swapped_tree, sequences, substitution_rate);

        if likelihood > best_likelihood {
            *tree = swapped_tree;
            best_likelihood = likelihood;
        }
    }
}

// 简单的子树交换函数
fn swap_subtrees(mut tree: TreeNode) -> TreeNode {
    if tree.left.is_some() && tree.right.is_some() {
        let temp = tree.left.clone();
        tree.left = tree.right.clone();
        tree.right = temp;
    }
    tree
}

// 计算节点的似然值
fn calculate_likelihood(node: &TreeNode, sequences: &HashMap<String, String>, substitution_rate: f64) -> f64 {
    // 如果是叶节点，直接返回
    if node.left.is_none() && node.right.is_none() {
        return 1.0;
    }

    // 获取左右子节点的似然值
    let left_likelihood = calculate_likelihood(node.left.as_ref().unwrap(), sequences, substitution_rate);
    let right_likelihood = calculate_likelihood(node.right.as_ref().unwrap(), sequences, substitution_rate);

    // 根据替换模型计算父节点的似然值
    let distance = node.distance;
    let transition_prob = f64::exp(-substitution_rate * distance);

    left_likelihood * right_likelihood * transition_prob
}

// 将树转换为 Newick 格式
fn tree_to_newick(node: &TreeNode) -> String {
    if node.left.is_none() && node.right.is_none() {
        return node.name.clone();
    }
    let left_str = tree_to_newick(node.left.as_ref().unwrap());
    let right_str = tree_to_newick(node.right.as_ref().unwrap());
    format!("({}:{},{}:{})", left_str, node.distance, right_str, node.distance)
}