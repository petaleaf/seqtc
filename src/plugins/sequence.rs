use crate::core::Plugin;
use clap::{ArgMatches, Command};
// use std::collections::HashMap;
pub struct SequencePlugin;

impl Plugin for SequencePlugin {
    fn name(&self) -> &'static str {
        "sequence"
    }

    // 定义该模块的使用参数
    fn command(&self) -> Command {
        Command::new("sequence")
            .about("获取DNA序列的反向互补序列")
            .arg(
                clap::arg!(-r --reverse "计算反向互补序列")
                    .required(false),
            )
            .arg(
                clap::arg!(-g --gc "计算GC含量")
                    .required(false),
                    
            )
            .group(
                clap::ArgGroup::new("structure") // 参数组名
                    .args(&["poly", "repeat"]) // 将 -p 和 -t 作为互斥参数
                    .required(false), // 可选互斥参数
            )
            // 分别定义 -p 和 -t 参数
            .arg(
                clap::arg!(-p --poly "检测poly结构")
                    .required(false), // 可选的poly结构参数
            )
            .arg(
                clap::arg!(-t --repeat "检测连续重复结构")
                    .required(false), // 可选的检测连续重复参数
            )
            .arg(
                clap::arg!(-b --base <BASE> "指定要检测的碱基")
                    .required(false), // 可选的指定碱基参数
            )
            .arg(
                clap::arg!(-m --min <MIN> "检测poly结构的最小重复次数/检测重复结构的最小重复次数")
                    .required(false)
                    .default_value("4"), // 默认最小重复次数为4
            )
            .arg(clap::arg!(<SEQUENCE> "DNA 序列"))

    }

    fn run(&self, matches: &ArgMatches) {
        let sequence = matches.get_one::<String>("SEQUENCE").expect("SEQUENCE is required");

        // 如果指定了 -r 或 --reverse 参数，则计算反向互补序列
        if matches.get_flag("reverse") {
            let reverse_complement = get_reverse_complement(sequence);
            println!("The reverse complement of {} is {}", sequence, reverse_complement);
        }

        // 如果指定了 -g 或 --gc 参数，则计算GC含量
        if matches.get_flag("gc") {
            let gc_content = calculate_gc_content(sequence);
            println!("The GC content of {} is {:.2}%", sequence, gc_content * 100.0);
        }
        // 检测poly结构
        if matches.get_flag("poly") {
            // 将输入的String解析为非负数
            let min_repeats: usize = matches.get_one::<String>("min")
                .unwrap()
                .parse()
                .expect("Invalid minimum repeat number");

            if let Some(base) = matches.get_one::<String>("base") {
                // 如果用户指定了要检测的碱基
                let poly_result = detect_poly_structure(sequence, base.chars().next().unwrap(), min_repeats);
                if poly_result.is_empty() {
                    println!("No poly structure detected for base {} with minimum {} repeats.", base, min_repeats);
                } else {
                    for (base, repeats, start_pos) in poly_result {
                        println!("Poly structure: Base {}, Repeats {}, Start Position {}", base, repeats, start_pos);
                    }
                }
            } else {
                // 没有指定碱基，检测所有可能的 poly 结构
                let all_poly_results = detect_all_poly_structures(sequence, min_repeats);
                if all_poly_results.is_empty() {
                    println!("No poly structures detected with minimum {} repeats.", min_repeats);
                } else {
                    for (base, repeats, start_pos) in all_poly_results {
                        println!("Poly structure: Base {}, Repeats {}, Start Position {}", base, repeats, start_pos);
                    }
                }
            }
        }

        // 检测连续重复结构
        if matches.get_flag("repeat") {
            let min_repeats: usize = matches.get_one::<String>("min")
            .unwrap()
            .parse()
            .expect("Invalid minimum repeat number"); // 获取最小重复数
            let repeats = detect_repeats(sequence, min_repeats);
            if repeats.is_empty() {
                println!("No repeating structures found.");
            } else {
                for (base, count, start_pos) in repeats {
                    println!(
                        "Repeat structure detected: Base {}, repeats {} times, starts at position {}",
                        base, count, start_pos
                    );
                }
            }
        }

        // 如果用户没有指定任何参数，则提示用户选择功能
        if !matches.get_flag("reverse") && !matches.get_flag("gc") {
            println!("Please specify either --reverse (-r) or --gc (-g) to perform an operation.");
        }
    }
}


fn get_reverse_complement(dna: &str) -> String {
    // 将序列转化为迭代器，这里可以接收实现了Deref<Target=str>的类型（&String也可以）
    dna.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c, // 处理非标准碱基
        })
        .collect()
}

fn calculate_gc_content(dna: &str) -> f64 {
    // 计算序列中GC的数
    let gc_count = dna.chars().filter(|&c| c == 'G' || c == 'C').count();
    gc_count as f64 / dna.len() as f64
}

fn detect_poly_structure(dna: &str, base: char, min_repeats: usize) -> Vec<(char, usize, usize)> {
    let mut result = Vec::new();
    let mut current_streak = 0;
    let mut start_pos = 0;

    for (i, c) in dna.chars().enumerate() {
        if c == base {
            if current_streak == 0 {
                start_pos = i; // 记录起始位置
            }
            current_streak += 1;
            if current_streak >= min_repeats {
                result.push((base, current_streak, start_pos));
            }
        } else {
            current_streak = 0; // 重新计数
        }
    }
    result
}

// 检测所有碱基的 poly 结构，返回符合条件的结构
fn detect_all_poly_structures(dna: &str, min_repeats: usize) -> Vec<(char, usize, usize)> {
    let mut results = Vec::new();
    let bases = vec!['A', 'T', 'C', 'G'];

    for base in bases {
        let base_results = detect_poly_structure(dna, base, min_repeats);
        results.extend(base_results);
    }
    results
}

fn detect_repeats(dna: &str, min_repeats: usize) -> Vec<(String, usize, usize)> {
    let mut repeats = Vec::new();
    let len = dna.len();

    // 遍历所有可能的子序列长度
    for sub_len in 2..=(len / min_repeats) {
        let mut i = 0;

        while i + sub_len <= len {
            let sub_seq = &dna[i..i + sub_len];
            let mut count = 1;
            let mut j = i + sub_len;

            // 检测子序列的连续重复
            while j + sub_len <= len && &dna[j..j + sub_len] == sub_seq {
                count += 1;
                j += sub_len;
            }

            if count >= min_repeats {
                repeats.push((sub_seq.to_string(), count, i));
            }

            // 继续查找下一个可能的起始位置
            i += 1;
        }
    }

    repeats
}