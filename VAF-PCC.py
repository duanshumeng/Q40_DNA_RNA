import pandas as pd
import numpy as np

# 生成有序编号列
column1 = list(range(1, 11))

# 生成随机的0到1之间的浮点数，5列10行
data = np.random.rand(10, 4)

# 创建数据框
df = pd.DataFrame(data, columns=['Column2', 'Column3', 'Column4', 'Column5'])
df.insert(0, 'Column1', column1) # 将有序编号列插入到第一列

# 计算皮尔森相关系数矩阵
correlation_matrix = df.iloc[:, 1:].corr()

# 提取第二列至第五列两两之间的皮尔森相关系数
correlation_pairs = correlation_matrix.unstack().sort_values(ascending=False)
correlation_pairs = correlation_pairs[correlation_pairs != 1] # 去除自身相关系数

# 打印数据框
print(df)

# 打印两两之间的皮尔森相关系数
print("Pearson Correlation Coefficients:")
print(correlation_pairs)