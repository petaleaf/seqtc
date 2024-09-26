use crate::core::Plugin;
use clap::{ArgMatches, Command};

pub struct AlignmentPlugin;

impl Plugin for AlignmentPlugin {
    fn name(&self) -> &'static str {
        "alignment"
    }

    fn command(&self) -> Command {
        Command::new("alignment")
            .about("对两个DNA序列进行比对")
            .arg(clap::arg!(<SEQ1> "第一个DNA序列"))
            .arg(clap::arg!(<SEQ2> "第二个DNA序列"))
    }

    fn run(&self, matches: &ArgMatches) {
        // 使用 get_one 获取参数值
        let seq1 = matches.get_one::<String>("SEQ1").expect("SEQ1 is required");
        let seq2 = matches.get_one::<String>("SEQ2").expect("SEQ2 is required");
        
        let (nw_score, aligned_a, aligned_b) = needleman_wunsch(seq1, seq2, 1, -1);
        println!("Needleman-Wunsch Score: {}", nw_score);
        println!("Aligned:\n{}\n{}", aligned_a, aligned_b);

        let (sw_score, aligned_a, aligned_b) = smith_waterman(seq1, seq2, 1, -1);
        println!("Smith-Waterman Score: {}", sw_score);
        println!("Aligned:\n{}\n{}", aligned_a, aligned_b);
    }
}

// Needleman-Wunsch 算法实现(全局比对)
fn needleman_wunsch(seq_a: &str, seq_b: &str, match_score: i32, gap_penalty: i32) -> (i32, String, String) {
    let m = seq_a.len();
    let n = seq_b.len();
    
    let mut score_matrix = vec![vec![0; n + 1]; m + 1];

    for i in 1..=m {
        score_matrix[i][0] = score_matrix[i - 1][0] + gap_penalty;
    }
    for j in 1..=n {
        score_matrix[0][j] = score_matrix[0][j - 1] + gap_penalty;
    }

    for i in 1..=m {
        for j in 1..=n {
            let match_value = if seq_a.chars().nth(i - 1) == seq_b.chars().nth(j - 1) {
                match_score
            } else {
                -match_score
            };

            score_matrix[i][j] = *[
                score_matrix[i - 1][j - 1] + match_value,
                score_matrix[i - 1][j] + gap_penalty,
                score_matrix[i][j - 1] + gap_penalty,
            ]
            .iter()
            .max()
            .unwrap();
        }
    }

    let mut aligned_a = String::new();
    let mut aligned_b = String::new();
    let mut i = m;
    let mut j = n;

    while i > 0 || j > 0 {
        if i > 0 && j > 0 && (score_matrix[i][j] == score_matrix[i - 1][j - 1] + if seq_a.chars().nth(i - 1) == seq_b.chars().nth(j - 1) { match_score } else { -match_score }) {
            aligned_a.push(seq_a.chars().nth(i - 1).unwrap());
            aligned_b.push(seq_b.chars().nth(j - 1).unwrap());
            i -= 1;
            j -= 1;
        } else if i > 0 && score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty {
            aligned_a.push(seq_a.chars().nth(i - 1).unwrap());
            aligned_b.push('-');
            i -= 1;
        } else {
            aligned_a.push('-');
            aligned_b.push(seq_b.chars().nth(j - 1).unwrap());
            j -= 1;
        }
    }

    aligned_a = aligned_a.chars().rev().collect();
    aligned_b = aligned_b.chars().rev().collect();

    (score_matrix[m][n], aligned_a, aligned_b)
}

// Smith-Waterman 算法实现（局部比对）
fn smith_waterman(seq_a: &str, seq_b: &str, match_score: i32, gap_penalty: i32) -> (i32, String, String) {
    let m = seq_a.len();
    let n = seq_b.len();
    
    let mut score_matrix = vec![vec![0; n + 1]; m + 1];
    let mut max_score = 0;
    let (mut max_i, mut max_j) = (0, 0);

    for i in 1..=m {
        for j in 1..=n {
            let match_value = if seq_a.chars().nth(i - 1) == seq_b.chars().nth(j - 1) {
                match_score
            } else {
                -match_score
            };

            score_matrix[i][j] = *[
                0,
                score_matrix[i - 1][j - 1] + match_value,
                score_matrix[i - 1][j] + gap_penalty,
                score_matrix[i][j - 1] + gap_penalty,
            ]
            .iter()
            .max()
            .unwrap();

            if score_matrix[i][j] > max_score {
                max_score = score_matrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    let mut aligned_a = String::new();
    let mut aligned_b = String::new();
    let mut i = max_i;
    let mut j = max_j;

    while i > 0 && j > 0 && score_matrix[i][j] > 0 {
        if score_matrix[i][j] == score_matrix[i - 1][j - 1] + if seq_a.chars().nth(i - 1) == seq_b.chars().nth(j - 1) { match_score } else { -match_score } {
            aligned_a.push(seq_a.chars().nth(i - 1).unwrap());
            aligned_b.push(seq_b.chars().nth(j - 1).unwrap());
            i -= 1;
            j -= 1;
        } else if score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty {
            aligned_a.push(seq_a.chars().nth(i - 1).unwrap());
            aligned_b.push('-');
            i -= 1;
        } else {
            aligned_a.push('-');
            aligned_b.push(seq_b.chars().nth(j - 1).unwrap());
            j -= 1;
        }
    }

    aligned_a = aligned_a.chars().rev().collect();
    aligned_b = aligned_b.chars().rev().collect();

    (max_score, aligned_a, aligned_b)
}