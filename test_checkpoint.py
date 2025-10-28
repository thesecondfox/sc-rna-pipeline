#!/usr/bin/env python3
"""
断点续传功能测试脚本

演示如何使用CheckpointManager进行断点续传
"""

import os
import json
import time
import random
from sc_pipeline_enhanced import CheckpointManager


def simulate_step(step_name, success_rate=0.8):
    """
    模拟流程步骤执行

    参数:
    -------
    step_name: str
        步骤名称
    success_rate: float
        成功概率（用于模拟随机失败）
    """
    print(f"\n{'='*60}")
    print(f"执行: {step_name}")
    print(f"{'='*60}")

    # 模拟工作负载
    for i in range(5):
        time.sleep(0.5)
        print(f"  处理中... {(i+1)*20}%")

    # 随机失败
    if random.random() > success_rate:
        raise Exception(f"{step_name} 模拟失败（用于测试断点续传）")

    print(f"✅ {step_name} 完成")


def test_checkpoint_basic():
    """测试基本的断点保存和恢复"""
    print("\n" + "="*80)
    print("测试1: 基本断点保存和恢复")
    print("="*80)

    output_dir = "./test_checkpoint_output"
    os.makedirs(output_dir, exist_ok=True)

    # 创建断点管理器
    ckpt = CheckpointManager(output_dir, auto_resume=True)

    # 模拟流程步骤
    steps = [
        ('step1_read', '数据读取'),
        ('step2_qc', '质量控制'),
        ('step3_integrate', '数据整合'),
        ('step4_annotate', '细胞注释')
    ]

    for step_id, step_name in steps:
        if not ckpt.is_completed(step_id):
            try:
                simulate_step(step_name, success_rate=1.0)  # 100%成功

                # 保存断点
                ckpt.save_checkpoint(
                    step_id,
                    file=f"{output_dir}/{step_id}.txt",
                    timestamp=time.time()
                )

                # 创建模拟文件
                with open(f"{output_dir}/{step_id}.txt", 'w') as f:
                    f.write(f"{step_name} 完成\n")

            except Exception as e:
                ckpt.mark_step_failed(step_id, e)
                print(f"\n❌ 步骤失败: {e}")
                break
        else:
            print(f"\n✓ {step_name} 已完成，跳过")

    # 显示进度
    print("\n")
    ckpt.print_status()

    print("\n✅ 测试1完成\n")


def test_checkpoint_with_failure():
    """测试失败后的断点恢复"""
    print("\n" + "="*80)
    print("测试2: 失败后断点恢复")
    print("="*80)

    output_dir = "./test_checkpoint_failure"
    os.makedirs(output_dir, exist_ok=True)

    # 清空之前的断点
    checkpoint_file = os.path.join(output_dir, '.checkpoint.json')
    if os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)

    steps = [
        ('step1_read', '数据读取'),
        ('step2_qc', '质量控制'),
        ('step3_integrate', '数据整合'),
        ('step4_annotate', '细胞注释')
    ]

    # 第一次运行：在步骤3失败
    print("\n>>> 第一次运行（步骤3会失败）")
    ckpt = CheckpointManager(output_dir, auto_resume=True)

    for idx, (step_id, step_name) in enumerate(steps):
        if not ckpt.is_completed(step_id):
            try:
                # 步骤3设置为失败
                success_rate = 0.0 if idx == 2 else 1.0
                simulate_step(step_name, success_rate=success_rate)

                ckpt.save_checkpoint(
                    step_id,
                    file=f"{output_dir}/{step_id}.txt"
                )

                with open(f"{output_dir}/{step_id}.txt", 'w') as f:
                    f.write(f"{step_name} 完成\n")

            except Exception as e:
                ckpt.mark_step_failed(step_id, e)
                print(f"\n❌ 步骤失败: {e}")
                print(f"💡 可以修复问题后继续运行\n")
                break
        else:
            print(f"\n✓ {step_name} 已完成，跳过")

    # 显示当前状态
    ckpt.print_status()

    # 第二次运行：从断点恢复
    print("\n>>> 第二次运行（从断点恢复）")
    ckpt = CheckpointManager(output_dir, auto_resume=True)

    # 显示恢复状态
    ckpt.print_status()

    # 继续执行剩余步骤
    for step_id, step_name in steps:
        if not ckpt.is_completed(step_id):
            try:
                simulate_step(step_name, success_rate=1.0)  # 这次都成功

                ckpt.save_checkpoint(
                    step_id,
                    file=f"{output_dir}/{step_id}.txt"
                )

                with open(f"{output_dir}/{step_id}.txt", 'w') as f:
                    f.write(f"{step_name} 完成\n")

            except Exception as e:
                ckpt.mark_step_failed(step_id, e)
                print(f"\n❌ 步骤失败: {e}")
                break
        else:
            print(f"\n✓ {step_name} 已完成，跳过")

    # 显示最终状态
    print("\n")
    ckpt.print_status()

    print("\n✅ 测试2完成\n")


def test_checkpoint_validation():
    """测试断点验证功能"""
    print("\n" + "="*80)
    print("测试3: 断点完整性验证")
    print("="*80)

    output_dir = "./test_checkpoint_validation"
    os.makedirs(output_dir, exist_ok=True)

    # 创建初始断点
    ckpt = CheckpointManager(output_dir, auto_resume=True)

    # 保存一个断点
    test_file = f"{output_dir}/test_data.txt"
    with open(test_file, 'w') as f:
        f.write("Test data\n")

    ckpt.save_checkpoint(
        'step1_read',
        file=test_file,
        n_cells=1000
    )

    print("\n✓ 创建了有效的断点")
    ckpt.print_status()

    # 模拟文件损坏（删除文件）
    print("\n>>> 模拟文件损坏（删除数据文件）")
    os.remove(test_file)

    # 重新加载断点管理器
    ckpt = CheckpointManager(output_dir, auto_resume=True)
    print("\n⚠️  检测到文件缺失")
    ckpt.print_status()

    print("\n✅ 测试3完成\n")


def test_checkpoint_reset():
    """测试断点重置功能"""
    print("\n" + "="*80)
    print("测试4: 断点重置")
    print("="*80)

    output_dir = "./test_checkpoint_reset"
    os.makedirs(output_dir, exist_ok=True)

    # 创建一些断点
    ckpt = CheckpointManager(output_dir, auto_resume=True)

    for i in range(1, 4):
        step_id = f"step{i}_test"
        ckpt.save_checkpoint(
            step_id,
            file=f"{output_dir}/{step_id}.txt",
            completed=True
        )

    print("\n已创建多个断点:")
    ckpt.print_status()

    # 重置断点
    print("\n>>> 重置所有断点")
    ckpt.reset(confirm=True)

    # 验证重置
    ckpt = CheckpointManager(output_dir, auto_resume=True)
    print("\n重置后:")
    ckpt.print_status()

    print("\n✅ 测试4完成\n")


def main():
    """运行所有测试"""
    print("\n" + "="*80)
    print("CheckpointManager 功能测试")
    print("="*80)

    tests = [
        ("基本断点保存和恢复", test_checkpoint_basic),
        ("失败后断点恢复", test_checkpoint_with_failure),
        ("断点完整性验证", test_checkpoint_validation),
        ("断点重置", test_checkpoint_reset)
    ]

    print("\n可用的测试:")
    for idx, (name, _) in enumerate(tests, 1):
        print(f"  {idx}. {name}")

    print(f"  0. 运行所有测试")

    choice = input("\n选择测试 (0-4): ").strip()

    if choice == '0':
        for name, test_func in tests:
            print(f"\n\n{'#'*80}")
            print(f"# {name}")
            print(f"{'#'*80}")
            test_func()
    elif choice in ['1', '2', '3', '4']:
        idx = int(choice) - 1
        name, test_func = tests[idx]
        print(f"\n\n{'#'*80}")
        print(f"# {name}")
        print(f"{'#'*80}")
        test_func()
    else:
        print("无效的选择")
        return

    print("\n" + "="*80)
    print("所有测试完成！")
    print("="*80)
    print("\n清理测试文件:")
    print("  rm -rf ./test_checkpoint_*")
    print()


if __name__ == "__main__":
    main()
