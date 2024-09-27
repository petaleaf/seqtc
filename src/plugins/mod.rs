// Module: plugins
// 声明插件
pub mod sequence;
pub mod alignment;

// 引入核心
use crate::core::Core;
use crate::plugins::sequence::SequencePlugin;
use crate::plugins::alignment::AlignmentPlugin;

// 安装插件
pub fn register_plugins(core: &mut Core) {
    core.register_plugin(Box::new(SequencePlugin));
    core.register_plugin(Box::new(AlignmentPlugin));
}
