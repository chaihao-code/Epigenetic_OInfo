# Reference identifiability report (Step 5)

模拟重复/条件: 50; 噪声 sigma ∈ {0.02, 0.05}; 面板 A/B × 模型 A/B

## Panel B / Model B
- bias@5% immune, sigma=0.05: -0.0017 (sd 0.0026); 检出率(>5%): 0.28
- bias@0% immune (假阳性): 0.0007
## Panel B / Model A
- bias@5% immune, sigma=0.05: 0.0016 (sd 0.0031); 检出率(>5%): 0.67
- bias@0% immune (假阳性): 0.0013
## Panel A / Model B
- bias@5% immune, sigma=0.05: -0.0016 (sd 0.0046); 检出率(>5%): 0.43
- bias@0% immune (假阳性): 0.0019
## Panel A / Model A
- bias@5% immune, sigma=0.05: 0.0003 (sd 0.0035); 检出率(>5%): 0.53
- bias@0% immune (假阳性): 0.0022

## 备注
- 混合在 beta 层面线性近似, 非 read-level; 已按要求单独说明。
- T/B 细分在 marker 数少时可能不稳定, 关注 Model A (lymphoid 合并)。