![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739195029350-a612fabe-bc42-462d-b401-36d6ea0fba7f.png)

此处使用python版本的scRNA-seq处理工具scanpy，而不是R版本的seurat，因为seurat包安装繁杂

一，准备工作

1，python库的安装：  
新建1个环境sc-python

```python
mamba create -n sc-python -c conda-forge -y scanpy python-igraph leidenalg python=3.12
conda activate sc-python
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739178647732-7f44a3ea-f2ba-4f9a-9a08-019e98e4a635.png)

2，单细胞相关分析参考：biomamba资料  
[https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247492295&idx=1&sn=d4e81588c0ac2906849c8bc44d079209&chksm=9b3c9b97ac4b1281ab462f866538f096ba2473a3d1426a27d8d97ccc636bfd081e1ba8934dfd&scene=21#wechat_redirect](https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247492295&idx=1&sn=d4e81588c0ac2906849c8bc44d079209&chksm=9b3c9b97ac4b1281ab462f866538f096ba2473a3d1426a27d8d97ccc636bfd081e1ba8934dfd&scene=21#wechat_redirect)

二，AnnData数据结构理解

AnnData被设计用于矩阵样数据的处理，可以方便的获得矩阵行与列的索引。例如在scRNA-seq的数据中，每一行的数据对应为一个细胞表达的所有基因数据，每一列的数据对应为一个基因在所有细胞中的表达数据(这点与R中不同，Seurat对象的矩阵正好是AnnData的转置版)。不仅是表达矩阵，每个细胞（如样本来源、细胞类型等）及每个基因（如别名、基因ID等）还会有自己的注释信息。另外，考虑到单细胞矩阵的稀疏性、以及数据结构需要有用户友好型的特性，AnnData显得再合适不过。AnnData的基本结构可以看这个示意图：

![](https://cdn.nlark.com/yuque/0/2025/jpeg/33753661/1739179421238-89fcb972-48f9-42dd-af49-90939a37ea8d.jpeg)

AnnData的数据结构是：行是样本（cell），列是gene

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739194880416-7477f3ef-c0b8-4ab9-8567-a125449e7216.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739194888353-cf904ef2-ecf7-4e52-878d-d001887d3e59.png)

1，构建模拟Anndata对象并查看

```python
# 导包
import numpy as np
import pandas as pd
import anndata as ad #导入 anndata 库，并将其命名为 ad。anndata 是一个用于处理单细胞基因表达数据的库，提供了 AnnData 对象，用于存储和操作高维数据
from scipy.sparse import csr_matrix #从 scipy.sparse 模块中导入 csr_matrix 类。csr_matrix 是一种压缩稀疏行矩阵格式，用于高效存储和操作稀疏矩阵
print(ad.__version__)
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739184784485-300e09c9-3cf5-4f4f-a95b-6a777d921ec2.png)

```python
# 模拟一个矩阵构建AnnData对象
counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32) #将生成的随机矩阵转换为 scipy.sparse 模块中的压缩稀疏行矩阵（CSR 格式），并指定数据类型为 float32
temp_adata = ad.AnnData(counts) 
print(f"这个AnnData对象包括{temp_adata.shape[0]}行，与{temp_adata.shape[1]}列") #使用 anndata 库的 AnnData 类创建一个 AnnData 对象，并将稀疏矩阵 counts 作为数据存储在其中
# 稀疏矩阵存放在temp_adata.X中
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739185620106-4795d84e-bf03-4db2-bea6-22cf889c3b62.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739184864709-ab5b92d2-f54c-4586-bdb9-4ccac56a687f.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739184920748-0da3d06e-af35-4b94-84fa-a016af83fe7c.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739184817993-dd286018-7f83-42e4-bd38-568f19e9a73c.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739185053163-1f329ea8-88cd-4ca1-883c-8a07a0854e9f.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739185034012-49cba067-7ee0-41c5-86bc-cd003031f2ab.png)

主要的构建函数

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739185221446-17e25518-e929-4d8e-ac70-f09ea2200f2d.png)

查看这里的counts，其实就是csr的稀疏矩阵格式，符合要求

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739185404552-06983c18-7db6-4da4-a7ab-0263b902414d.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739185415456-621ca2ab-a302-43f8-a939-75ad30f85f15.png)



对于构建生成的这个注释之后的数据对象，进行仔细分析查看：

数据对象

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739185689451-43e5eb7a-c469-416b-afdb-53d92c2fdefe.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186357954-1164aef3-2690-4897-a1f8-61ccc53a5283.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739194903829-0c693fc6-7ab1-45ce-b2eb-754fa0ab4843.png)

数据类型

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739185736916-fc2ca8c3-880a-421e-accf-a84b93eeb264.png)

尺寸形状

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739185784603-f13b5e1c-fcd8-4801-bdcb-916b7e491057.png)



```python
蓝色盒子：

该图标通常表示一个标准的类或对象。它表示对象具有相关的方法和属性，并可以通过其定义的函数或属性进行交互。这一视觉提示暗示对象内部存在一个明确定义的结构。
紫色盒子：

通常表示特定类型的对象或命名空间，特别是模块或包。这个符号暗示对象是一个集合，可能包含多个函数或其他子对象，为代码中的层级关系提供了视觉指示。
白色方框：

通常表示一个普通的函数或变量。该图标通常显示基本信息，例如函数签名或类型提示，表明这是一个简单的对象，没有特殊的内部结构。它为开发者提供了基础层次的信息。
扳手：

常用于表示可以对对象执行的一系列操作，如重构选项或编辑功能。此图标可能表明重命名、导航至对象定义或应用其他修改等操作，从而促进代码维护。
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739185857027-ef09b584-6549-4ab0-bdd9-49d1bb0195b9.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186381220-03db2f82-34e0-4178-bee3-1b7cddeca5e2.png)

X存放的就是稀疏矩阵本身

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186472029-454b63f6-fd1e-4b9a-a749-7d5a7b0ce7b4.png)

行

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186616480-d2315aaa-002d-4cfd-b0e0-f40e2f50f6b8.png)

行数

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186896886-5ba6615e-0a49-4307-a9a0-fe70fb3a02ac.png)

行名

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187192110-3e0a3bad-a876-4dde-86a7-165f2e4b36a4.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187259164-476e4775-8f52-483a-875c-0b0e161f576f.png)

列

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186646318-0bafb410-519f-4455-b201-27bcd49cf4ac.png)

列数

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187069043-1d1eef39-b7e3-4de7-b2a3-d6c48aa4edb4.png)

列名

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187211063-5ffceef1-ad83-426a-bab4-f31d072201bf.png)



后面都是空的slot

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186663855-e234e8e4-1507-42a6-9a17-0f1c7cd8a09e.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186683288-b4a7eb7b-6954-4f32-a246-f23057ede56b.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186704505-89fe6920-036a-4c7a-a3c8-f5c86977402d.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186725637-8522b348-62e5-4be2-98b7-56bc04f4281a.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739186764212-8e88749b-2098-468f-9fa6-b6f90a8184e2.png)



AnnData索引操作可以帮助我们取得对应的AnnData子集，这样我们就可以选取感兴趣的细胞或基因参与下游的计算，简化流程并节省计算资源。类似于pandas的DataFrame通过obs_names和var_names、布尔值、索引整数数字均可以完成对AnnData的取子集操作。（切片slice，bool/index）

我们的矩阵是刚才模拟出来的，可以用这样的方式添加上细胞名称与基因名称：

行使观测obs样本cell，列是变量var基因

```python
# 生成细胞名称并传递给obs_names
temp_adata.obs_names = [f"Cell_{i:d}" for i in range(temp_adata.n_obs)]
# 生成基因名称并传递给var_names
temp_adata.var_names = [f"Gene_{i:d}" for i in range(temp_adata.n_vars)]

print(f"前五个细胞名称为:\n{temp_adata.obs_names[:5]}")
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187512792-3e7c4759-ac3d-4f77-b53d-f90eda975f02.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187317954-57d7c450-0804-47ca-8d33-e1b29b0adc9d.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187330358-0a99f1e1-e2fb-4409-b993-15a175d35236.png)

修改之后

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187409999-968aa455-eea0-4439-b631-c0651d75b0af.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187419417-afdbac93-8e60-4d9f-8dd1-9968dc322dd0.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187432075-4bcabf9d-e7e8-407c-896a-203c53b3fc23.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187449580-56ce30cd-807b-4385-9ffc-95a5723cc536.png)

就可以通过行列切片方式取出1个子集：bool或index

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187648645-d9c46cda-9430-42b4-b367-7eee4d3562cd.png)

```python
# 那么子集就可以这么取:
temp_adata[["Cell_1", "Cell_10"], ["Gene_5", "Gene_1900"]]
# 这样就得到了一个包含"两个细胞"、"两个基因"的表达矩阵
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187669623-0d84ea07-19f3-42b8-942f-169247a833e6.png)

2，添加注释信息

其实前面已经观察到obs以及var都是一个slot对象，两个都是pandas.Dataframe对象，熟悉pandas语法即可

obs是行对象数据框，var是列对象数据框

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739187841171-146ea003-22ac-435d-a1a1-e431ecc2f630.png)

adata.obs与adata.var均是DataFrame，那么就可以很容易的加上对应的注释，例如加上细胞类型的注释:



```python
# 例如这里随机生成"B", "T", "Monocyte"三种细胞类型作为注释
ct = np.random.choice(["B", "T", "Monocyte"], size=(temp_adata.n_obs,))
print(f"ct的前10个元素为{ct[0:10]}")

# 加入adata.obs中
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739188734086-937a433b-e0b3-4a18-b8bb-4a1e65f6a2f1.png)

生成的是一维数组

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739188672159-f4a3a72c-16cc-4f26-b398-8050578a625d.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739188725761-9397f843-254e-4a79-a022-3fd84716b1ff.png)



从提供的列表list（cell名称注释中）随机选取

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739188223759-789dc89b-190f-4a68-a77e-f74a19d3989f.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739188333314-466cd216-d461-43ca-809b-cb6a1da23774.png)

然后这里的size提供int还是元组都可以

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739188576046-44f953dd-c69a-44d2-9d11-984867259065.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739188582250-f91ed725-590d-4d4c-845e-06d42360e14f.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739188588440-29c14071-55d2-4272-964a-d9142d9c8ff4.png)

前面我们知道行obs是一个数据框，其实这个数据框没有列column，即var

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739188874361-2653a24f-a452-4f20-8d58-d1439fe8b6f7.png)

可以新增一列，也就是注释



```python
temp_adata.obs["cell_type"] = pd.Categorical(ct)  # Categorical的执行效率较高

# 可以看到adata.obs已经包含了我们刚刚添加的注释
temp_adata.obs
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739189280971-d98e415f-3740-4324-bd2d-c6f6497262b4.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739189319837-015b3910-5e2f-41f0-b094-f9ada9ca40b3.png)



因为细胞cell很多，但是一般注释出来的类型很少，就固定那么几种，故cell类型其实就是分类变量，即R中的factor因子类型变量，其实就是相当于在R中做了一层factor（）函数类型转换处理

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739189185537-a0826037-2453-4913-a5ef-4bf472fc4288.png)

返回一维数组

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739189245120-34cbbad0-fd7f-4ec9-a866-0380406016f6.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739189344808-26c69c48-1e3b-4c39-a50d-64985ce74df7.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739189045818-ba563554-adff-42b7-b185-56f4d8968537.png)

其实上面就相当于是给数据框添加列变量，因为行index已经提供，自己再补充填充列var即可

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739188037793-e4899542-eadb-4cd9-b682-ff7923f2d1b5.png)



同样的，对行数据框做注释（添加新列）之后，也可以对列数据框做注释（添加新列）![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739189464007-89aa7378-7856-40bc-96f8-48ce3cb4b992.png)

同样处理，先生成对应的注释一维数组

```python
# 例如这里随机生成"Pathway_A", "Pathway_B", "Pathway_C"三种通路名作为细胞的注释
my_path = np.random.choice(["Pathway_A", "Pathway_B", "Pathway_C"],
                      size=(temp_adata.n_vars,))
print(f"ct的前10个元素为{my_path[0:10]}")
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739189512075-b811b12f-6364-427b-bec4-e358d7195b41.png)

同样是对列数据框添加列注释

```python
temp_adata.var["my_pathway"] = pd.Categorical(my_path)  # Categorical的执行效率较高

# 可以看到adata.var已经包含了我们刚刚添加的注释
temp_adata.var
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739189594034-7f437921-6d38-4b2d-813e-67a3cf69cb6b.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739189613870-028e1ea9-dfa2-422c-8230-88401ee75d01.png)

行对象数据框obs的注释cell_type以及列对象数据框var的注释my_pathway;

然后这个时候可以进一步使用行以及列的注释来提取子集：

比如说对于行索引，行是cell，其属性（var）目前只有1列，即cell类型，那我就通过行对象的这一个属性来提取或者说是范围限制特定的行，

即  数据【数据行的条件，：】用于提取符合某些行条件的子集

其中**数据行的条件 **为 **行对象数据框.列属性=某某条件**

```python
# 这时可以通过注释来选取数据的子集：
# 取出B cell的细胞数据
bcell_adata = temp_adata[temp_adata.obs.cell_type == "B",:]

# 取出仅包含Pathway_A的基因数据：
PA_adata = temp_adata[:,temp_adata.var.my_pathway == "Pathway_A"]
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739190272804-e78bdab9-a8f7-4fa8-9936-8072e60db474.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739190034732-a9051e79-63c3-4ade-a78a-f55a116a3401.png)  
![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739190114056-3327162f-a8f7-40eb-80ec-bc2e4cd073e0.png)

在该单细胞的稀疏矩阵中，属于B类型的只有40个细胞，总共是100个细胞；

属于通路A的只有662个gene，总共是2000个gene；



当然，在实际分析过程中，往往会有非常复杂且大量的注释，比如说对于行对象，也就是cell，

可能会有时序信息，也就是发育等信息；

病患来源信息，以及亚型信息，以及组织采集位点信息等

```python
# 构建一个略显"复杂"的注释DataFrame
obs_meta = pd.DataFrame({
        'time_yr': np.random.choice([0, 2, 4, 8], temp_adata.n_obs),
        'subject_id': np.random.choice(['subject 1', 'subject 2', 'subject 4', 'subject 8'], temp_adata.n_obs),
        'instrument_type': np.random.choice(['type a', 'type b'], temp_adata.n_obs),
        'site': np.random.choice(['site x', 'site y'], temp_adata.n_obs),
    },
    index=temp_adata.obs.index,   # 需要与AnnData的observations相一致
)
obs_meta.head()
```

创新数据框的时候需要提供行索引temp_adata.obs.index

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739190788850-4eb4384f-f60b-4331-92d7-346f72fcd3da.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739190807529-126f0348-922e-4c16-845a-6fae3176ca67.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739190816028-7168113d-7914-4c48-b668-84e496b43fe1.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739190824409-b866a01a-1233-4bc9-83e8-4c8c19ac02ec.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739190832591-4174a169-c9b1-4089-a476-cae29d63502f.png)

而这仅仅只是行的注释，即对cell的注释，实际过程中还会有对gene的注释

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739190880426-dc031b13-dd35-419b-ba65-3e050643f7bd.png)



现在再回到我们这里的定义：

我们在定义这个注释数据框的时候可以提供很多对象以及注释信息

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739190961005-b67caf00-c4e3-41c1-893e-80f92388d5e6.png)

比如说是提供稀疏矩阵，以及分别行数据框以及列数据框对象的注释信息

```python
# 将注释添加到AnnData对象中:
new_adata = ad.AnnData(temp_adata.X, obs=obs_meta, var=temp_adata.var)

new_adata
# 可以看到多出很多注释内容
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739191105088-ac222649-30d4-4f2a-acfa-d1852ebd6fe2.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739191111079-bf58d031-46b8-44bf-a353-215cf48a62e1.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739191123644-8402ce89-668b-4b56-ac9e-b69b308dd9b0.png)

3，注释矩阵管理

刚才我们添加的cell_type以及my_pathway包含的均是单个维度的注释信息（因为我们都是使用np.random.choice随机生成的一维矩阵，都是一维的注释信息），

而有些数据是以多维度的形式参与单细胞的数据处理。例如UMAP与tsne的降维矩阵，这些多维度的信息可以添加至.obsm或.varm中。需要注意的是，矩阵的长度需要与.n_obs和.n_vars相对应。

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739192107035-8db98911-6d08-4150-a72a-1cddaad21d86.png)

从上面我们可以看出，一维的注释信息我们可以从obs以及var中直接新增一列，

如果是多维的数据，我们需要在obsm以及varm中添加（m指代的是multi维度）；

然后同样的，我们将obs以及var当做一维的注释数据框对象，

obsm以及varm也可以直接用于当做多维的注释数据框对象，我们同样可以像新增列名一样去增加属性

例如这里我们模拟一个UMAP和gene_stuff的结果进行添加：

```python
temp_adata.obsm["X_umap"] = np.random.normal(0, 1, size=(temp_adata.n_obs, 2))
temp_adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(temp_adata.n_vars, 5))
temp_adata
# 可以看到X_umap与gene_stuff数据被添加
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739192630350-035c80a4-8293-4d2e-9b65-f769e4e65f8e.png)

其实就是行cell多了一个2维的注释信息，列gene多了一个5维的注释信息

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739192523244-30aaf3e5-f815-47d7-9f27-7308a0c41187.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739192589584-cee9a5d3-309e-4356-a5be-a7cb07722d60.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739192603331-11674c6b-2116-4c66-9b2e-924b47825fd2.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739192617776-a01d3aef-ec5c-4cea-bf84-243c387dbced.png)

其实和前面一样的道理，一维的数据我们直接对行obs以及列var对象新增一列，增加一维数组即可；

多为的数据，我们暂时没法直接在obs以及var对象中添加，所以需要在其他的参数对象中添加



4，非结构性的注释信息

这部分Unstructure metaData可以存放在.uns中，并且对数据格式、内容均无任何要求，可以存放一些提示信息：

同样是参数uns+新增1列var

```python
temp_adata.uns["random"] = [1, 2, 3]
temp_adata.uns
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739192825955-869989da-4472-4abb-9eb3-525dc140737d.png)



5，Layers

原始的矩阵数据可能会经计算生成一些新的数据，这时我们就可以将这些过程或结果数据存放在不同的Layers中，例如数据处理过程中产生的标准化数据log_transformed；

其实我个人觉得可以理解为channel也就是通道，因为我们的表达矩阵需要进行归一化等数据处理，所以相当于是不同层的矩阵数据叠加在一起，当然是不同的channel；

然后行或者是列的注释数据其实一般是不会变的，因为矩阵只是值发生了映射变换，实际上并没有变形shape，所以我们的行或者是列的注释信息还是对应的

```python
temp_adata.layers["log_transformed"] = np.log1p(temp_adata.X)
temp_adata
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739193806113-c1d8cf2a-0fe8-4549-a2e4-1d6096422aa4.png)

注意log1p是自然对数e，是自然natural

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739193198373-48d4791d-df67-41b7-b610-e647259d6337.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739193783602-192a2366-984a-48a1-ae3a-acacdae22f48.png)

先在layer中定义1个新矩阵，即先定义1个layer，再将layer转换为数据框dataframe



```python
temp_adata.to_df(layer="log_transformed")# Layer也可以直接转换为矩阵
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739194100674-5c9c776c-e665-4d67-827d-56ae1ae726a6.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739194046358-4c73dd79-43b2-450b-998f-1a47d4f07848.png)

将layer转换为df

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739194062773-d8a60571-19a1-4c43-9f7d-5cb6be4b60c6.png)

```python
# 也可以输出为csv文件：当然前面转换为数据框之后可以直接导出为csv文件
temp_adata.to_df(layer="log_transformed").to_csv("log_transformed.csv")
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739194511756-af26a0af-b160-4100-af85-7dd867327eda.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739194390329-a0aea63a-76f1-4df6-a448-e48e2ecc08b1.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739194452634-967fa0d1-b496-4492-ba70-8e409b5b1649.png)

6，AnnData存储与读取

AnnData可以很方便的以h5ad的形式保存在本地，这种格式的读写速度实际体验起来均比rds文件要快。

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739195003087-1741d6fc-dc7f-48d2-b974-d16d5519e7fb.png)

```python
#保存为h5ad格式
temp_adata.write('my_results.h5ad', compression="gzip")
```

实际上就是AnnData object对象的write属性而已

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739195408446-70ea1123-87ab-462e-8d95-33d54082cb89.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739195529669-5095ff65-4e56-42eb-8828-a4b226a10209.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739195564567-69aabecb-4140-415b-b0c9-5c752c16a4f9.png)



```python
adata = ad.read('my_results.h5ad', backed='r')
adata.isbacked #此时my_results.h5ad处于open状态，即被占用的状态
```

FutureWarning: `anndata.read` is deprecated, use `anndata.read_h5ad` instead. `ad.read` will be removed in mid 2024

根据警告信息，anndata.read 方法已被弃用，建议使用 anndata.read_h5ad 方法来读取 .h5ad 文件

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739195758173-2ccdbf90-4267-4213-8c92-731007c58a74.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739196202373-1355c53d-c8be-4669-8be8-aafd13aa84e4.png)

```python
# 从文件中读取 AnnData 对象，以只读模式打开
adata = ad.read_h5ad('my_results.h5ad', backed='r')

# 检查 AnnData 对象是否处于 backed 模式
print(adata.isbacked)  # 输出: True

# 关闭文件，释放资源
adata.file.close()
print(adata.isbacked)  # 输出: True
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739196262540-6f4ef31a-8ed2-4919-a57d-0b30eb8b62b6.png)

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739196370094-133e48ac-dfe0-41ae-994d-5cc48bb0892d.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739196525744-074300fb-bd66-44ee-9d95-c8356dd43542.png)



```python
adata = ad.read_h5ad('my_results.h5ad', backed='r')
adata.isbacked #此时my_results.h5ad处于open状态，即被占用的状态，输出为True

# 此时adata可以被正常使用，且多出路径变量
adata.filename

# 当完成分析时可以解除占用：
adata.file.close() 
# 这其实就是scanpy在分析相同的数据比Seurat占用更少内存的秘诀
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739196608332-27caaa89-001f-4f10-8111-4abcd40fa2d2.png)



7，AnnData的Views与copies

与pandas的DataFrame一样，AnnData也可以通过Views来创建一个对象的”影子”，这样既可以避免占据新的内存，又可以修改原有的对象

```python
# 在做View时候也可以使用切片：
temp_adata[:5, ['Gene_1', 'Gene_3']]


# 而copy的方式会占据额外的内存
adata_subset = temp_adata[:5, ['Gene_1', 'Gene_3']].copy()
adata_subset
```

![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739196950994-4c2c5819-69d6-4649-bd11-f07534e06d7a.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739196968977-78fad229-d180-4037-b375-cbe20268d5a9.png)![](https://cdn.nlark.com/yuque/0/2025/png/33753661/1739196977148-cf800c7b-1a35-45ec-91b2-11753a58a652.png)

