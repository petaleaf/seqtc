use crate::core::Plugin;
use clap::{ArgMatches, Command};

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