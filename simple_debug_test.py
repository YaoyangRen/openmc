#!/usr/bin/env python3
"""
简单的测试脚本来调试格林函数问题
"""

import os
import sys
import subprocess

def create_simple_test():
    """创建一个简单的测试配置"""
    
    # 创建测试目录
    if not os.path.exists("simple_test"):
        os.makedirs("simple_test")
    
    os.chdir("simple_test")
    
    # 创建简单的几何文件
    geometry_xml = """<?xml version="1.0"?>
<geometry>
  <cell id="1" material="1" region="-1" />
  <cell id="2" material="void" region="1" />
  
  <surface id="1" type="sphere" coeffs="0.0 0.0 0.0 5.0" />
</geometry>
"""
    
    # 创建材料文件
    materials_xml = """<?xml version="1.0"?>
<materials>
  <material id="1">
    <density value="1.0" units="g/cm3" />
    <nuclide name="H1" ao="2.0" />
    <nuclide name="O16" ao="1.0" />
  </material>
</materials>
"""
    
    # 创建设置文件
    settings_xml = """<?xml version="1.0"?>
<settings>
  <run_mode>fixed source</run_mode>
  <particles>100</particles>
  <batches>2</batches>
  
  <source>
    <space type="point" parameters="0.0 0.0 0.0" />
    <angle type="isotropic" />
    <energy type="discrete" parameters="1.0e6" />
  </source>
</settings>
"""
    
    # 写入文件
    with open("geometry.xml", "w") as f:
        f.write(geometry_xml)
    
    with open("materials.xml", "w") as f:
        f.write(materials_xml)
    
    with open("settings.xml", "w") as f:
        f.write(settings_xml)
    
    print("简单测试配置已创建")

def run_simple_test():
    """运行简单测试"""
    print("运行简单测试...")
    
    try:
        result = subprocess.run(["../build/bin/openmc"], 
                              capture_output=True, 
                              text=True, 
                              timeout=60)
        
        print("返回码:", result.returncode)
        print("标准输出:")
        print(result.stdout)
        
        if result.stderr:
            print("标准错误:")
            print(result.stderr)
            
        return result.returncode == 0
        
    except Exception as e:
        print(f"运行失败: {e}")
        return False

def main():
    """主函数"""
    print("=" * 50)
    print("简单测试程序")
    print("=" * 50)
    
    original_dir = os.getcwd()
    
    try:
        create_simple_test()
        success = run_simple_test()
        
        if success:
            print("测试成功完成")
        else:
            print("测试失败")
            
    except Exception as e:
        print(f"错误: {e}")
        
    finally:
        os.chdir(original_dir)

if __name__ == "__main__":
    main()