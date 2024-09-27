// 使用clap实现插件化的核心
use clap::{ArgMatches, Command};

pub trait Plugin {
    fn name(&self) -> &'static str;
    fn command(&self) -> Command; // 插件定义自己的命令
    fn run(&self, matches: &ArgMatches); // 插件执行逻辑
}

pub struct Core {
    plugins: Vec<Box<dyn Plugin>>,  // 插件列表
}

/*
Core 结构体包含一个插件列表，每个插件都是一个实现了 Plugin 特征的 Box<dyn Plugin>。
他包含三个功能：
1. 添加插件
2. 运行插件
3. 列出所有插件
*/
impl Core {
    pub fn new() -> Self {
        Core { plugins: Vec::new() }
    }
    pub fn list_plugins(&self) {
        for plugin in &self.plugins {
            println!("Registered plugin: {}", plugin.name());
        }
    }

    pub fn register_plugin(&mut self, plugin: Box<dyn Plugin>) {
        self.plugins.push(plugin);
    }

    // 接收main函数传入的参数，根据参数执行对应的插件
    pub fn build_cli(&self) -> Command {
        let mut cli = Command::new("bio")
            .about("A bioinformatics framework with plugin support")
            .subcommand(
                Command::new("list")
                    .about("列出所有已注册的插件")
            );
        
        // 将每个插件的命令添加为子命令
        for plugin in &self.plugins {
            cli = cli.subcommand(plugin.command());
        }

        cli
    }

    pub fn run(&self, matches: &ArgMatches) {
        if let Some((command, sub_matches)) = matches.subcommand() {

            if command == "list" {
                self.list_plugins();
                return;
            }



            for plugin in &self.plugins {
                if plugin.name() == command {
                    plugin.run(sub_matches);
                    return;
                }
            }
        }
        println!("Unknown command: {}", matches.subcommand_name().unwrap_or("none"));
    }


}
