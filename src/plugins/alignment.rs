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
        
        let result = align_sequences(seq1, seq2);
        println!("Alignment result between {} and {}: {}", seq1, seq2, result);
    }
}

fn align_sequences(seq1: &str, seq2: &str) -> String {
    if seq1 == seq2 {
        "完全匹配".to_string()
    } else {
        "不匹配".to_string()
    }
}
