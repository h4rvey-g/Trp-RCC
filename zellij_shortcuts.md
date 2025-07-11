# Zellij 快速入门指南 (与 Tmux 对比)

欢迎来到 Zellij！如果您已经熟悉 tmux，您会发现很多相似之处，但 Zellij 旨在提供更直观、更友好的用户体验。

## 核心概念：锁定模式

与 `tmux` 的 `prefix` (默认为 `Ctrl-b`) 不同, Zellij 使用 **锁定模式**。

*   **进入锁定模式**: 按 `Ctrl-p`。一旦进入锁定模式，您可以按后续键来执行命令 (例如, `d` 来分离会话)。
*   **退出锁定模式**: 按 `Enter` 或 `Esc`。

这就像一个临时的命令模式，而不是每次操作前都需要按前缀键。

## 基本操作

| 功能 | Zellij 快捷键 | Tmux 快捷键 | 备注 |
| :--- | :--- | :--- | :--- |
| **分离会话** | `Ctrl-p` -> `d` | `Ctrl-b` -> `d` | "d" 代表 "detach"。 |
| **列出/切换会话** | `zellij list-sessions` / `zellij attach <name>` | `tmux ls` / `tmux attach -t <name>` | 在 shell 中执行。 |
| **打开命令帮助** | `Ctrl-p` -> `?` | `Ctrl-b` -> `?` | |

## 窗格 (Panes) 管理

Zellij 的窗格管理非常直观。

| 功能 | Zellij 快捷键 | Tmux 快捷键 | 备注 |
| :--- | :--- | :--- | :--- |
| **新建窗格 (向下)** | `Ctrl-p` -> `n` (然后 `d`或`↓`) 或 `Alt-n` | `Ctrl-b` -> `"` | 在 Zellij 中，`n` 是 `new`。 |
| **新建窗格 (向右)** | `Ctrl-p` -> `n` (然后 `r`或`→`) 或 `Alt-r` | `Ctrl-b` -> `%` | `r` 是 `right`。 |
| **在窗格间移动** | `Alt` + `↑` `↓` `←` `→` | `Ctrl-b` + `↑` `↓` `←` `→` | Zellij 的方式更直接。 |
| **关闭当前窗格** | `Ctrl-p` -> `x` 或 `exit` | `Ctrl-b` -> `x` 或 `exit`| `x` 代表 "exit"。 |
| **最大化/恢复窗格** | `Ctrl-p` -> `z` | `Ctrl-b` -> `z` | `z` 代表 "zoom"。 |
| **浮动/嵌入窗格** | `Ctrl-p` -> `w` | (无默认) | `w` 代表 "toggle floating"。|

## 标签页 (Tabs) 管理

| 功能 | Zellij 快捷键 | Tmux 快捷键 (Windows) | 备注 |
| :--- | :--- | :--- | :--- |
| **新建标签页** | `Ctrl-t` | `Ctrl-b` -> `c` | `t` 代表 "tab"。 |
| **切换到下一个标签页**| `Ctrl-t` (或 `Alt-l`, `Alt-→`) | `Ctrl-b` -> `n` | `n` for "next"。|
| **切换到上一个标签页**| `Ctrl-t` (或 `Alt-h`, `Alt-←`) | `Ctrl-b` -> `p` | `p` for "previous"。 |
| **直接切换到标签页 N** | (无默认, 可配置) | `Ctrl-b` -> `N` (数字) | |
| **重命名标签页** | `Ctrl-p` -> `r` | `Ctrl-b` -> `,` | `r` 代表 "rename"。|
| **关闭当前标签页** | `Ctrl-p` -> `q` | `Ctrl-b` -> `&` | `q` 代表 "quit"。|

## 其他实用功能

| 功能 | Zellij 快捷键 | Tmux 快捷键 | 备注 |
| :--- | :--- | :--- | :--- |
| **进入滚动模式** | `Ctrl-s` | `Ctrl-b` -> `[` | `s` for "scroll"。 |
| **在滚动模式中移动**| `↑` `↓` / `PageUp` `PageDown` | `↑` `↓` / `PageUp` `PageDown` |  |
| **退出滚动模式** | `Enter` or `Esc` | `q` | |
| **调整窗格大小** | `Ctrl-r` | `Ctrl-b` + `Ctrl-` + `↑` `↓` `←` `→` | `r` for "resize"。 |
| **锁定界面** | `Ctrl-g` | (无默认) | `g` for "go" (away)。冻结界面，防止意外输入。 |

## 如何开始

1.  **启动一个新的 Zellij 会话**:
    ```bash
    zellij
    ```
2.  **启动一个带名字的会话**:
    ```bash
    zellij -s my-session-name
    ```
3.  **重新连接 (Attach) 到最后一个会话**:
    ```bash
    zellij attach
    ```
4.  **连接到指定名字的会话**:
    ```bash
    zellij a my-session-name
    ```

这个列表涵盖了大部分日常使用所需的功能。Zellij 的一个巨大优势是其配置和布局系统，允许您保存和加载复杂的窗格布局。祝您使用愉快！