// src/core.rs
// 插件特征，定义插件的基本行为
pub trait Plugin {
    fn name(&self) -> &'static str;
    fn run(&self, args: Vec<String>);
}

pub struct Core {
    // 在堆上储存实现了 Plugin 特征的插件
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

    // 添加插件
    pub fn register_plugin(&mut self, plugin: Box<dyn Plugin>) {
        self.plugins.push(plugin);
    }
    // 运行插件
    pub fn run_plugin(&self, name: &str, args: Vec<String>) {
        for plugin in &self.plugins {
            if plugin.name() == name {
                plugin.run(args);
                return;
            }
        }
        println!("Plugin not found: {}", name);
    }
    // 列出所有插件
    pub fn list_plugins(&self) {
        for plugin in &self.plugins {
            println!("Registered plugin: {}", plugin.name());
        }
    }
}
