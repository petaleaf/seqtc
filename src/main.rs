mod core;
mod plugins;

// use clap::Parser;
use core::Core;

fn main() {
    let mut core = Core::new();
    
    // 注册插件
    plugins::register_plugins(&mut core);

    // 解析命令行输入
    let cli = core.build_cli().get_matches();

    // 根据解析结果运行插件
    core.run(&cli);
}