# My
This is a forked repo from https://github.com/InteractiveComputerGraphics/SPlisHSPlasH


这是对Jan Bender原本repo的fork。原repo为：https://github.com/InteractiveComputerGraphics/SPlisHSPlasH

所有的增改都尽量保证开闭原则，即宁增加，勿修改。尽量保证不动原有代码。凡是要增加新功能，都要新建一个类/自由函数，而不破坏原有代码。即使要修改原代码，也建议通过复制粘贴原代码，新建一个类，并在新的文件中更改。这样的好处是尽量让修改是可退的，保证出现差错之后，我们还能恢复到可以运行的版本。

我们对源码的修改都尽量几种放在SPlisHSPlasH/My这个文件夹下面。我们自己新增的场景都放在 data/MyScenes下面


我们在此基础上进行了许多增改：
1. Diffusion: 按照温度进行扩散。请看 temperatureDiffusion和Coagulation。
2. Viscosity_Casson: 请看 Viscosity/Viscosity_Casson。根据Casson公式计算黏度，并且根据温度动态调整黏度。
3. Plasticity: 弹塑性，请看塑性冲击的例子plastic_strike.json。请看 My/Plasticity。
4. (有BUG)ShapeMatching: 刚体。请看RigidBody_Bunny.json的例子。
5. (未完成)NonNewton: 请看 NonNewton。非牛顿剪切变稀和剪切变稠。
6. MyPartioReader：(目前位于Utilities文件夹)粒子的导入器，可以直接导入Houdini的 .bhclassic 文件(.bhclassic是Houdini旧版（12以前）的bgeo文件）。 可以读入uv
7. 用户交互：现在可以用鼠标右键旋转视角，鼠标滚轮缩放视角，鼠标中键平移视角。更加符合人类直觉。
8. Interactive: 交互类。可以获取点击的鼠标世界位置（打印到屏幕上）。还可以控制刚体（通过control+WASDF控制第0号刚体）。
9. SurfaceParticles: 自由函数findSurfaceParticles。用于找到所有表面粒子（定义为邻居数量<10的粒子）
10.  Exporter_xyz: 导出粒子信息为xyz格式（就是每行输出ascii）。目前位于Simulator/exporter中
11.  Exporter_MyPartio: 导出的bgeo格式可以直接被Houdini读取。修复了原有partio导出器的BUG。目前位于Simulator/exporter中。
12. (未完成)MyTimeStep: 留空，以后将DFSPH搬到这里修改。

## 非牛顿
非牛顿目前实现方式是通过NonNewton类，根据不同的子模型计算出粘度，然后将粘度付给不同的粘性子模型。这就需要对原有的粘性子模型进行非常轻微但侵入式的改动。基本上是将原本的粘度（一个标量）改为向量，让每个粒子的粘度都不相同。


## 其他

**导出自定义属性的方式**

在json中particleAttributes中给出自定义场的名称。该名称应该和构造函数中的model->addField的名称一致。

勾选Houdini Exporter(MyPartio)或在json中给定enableMyPartioExport为true。

**模拟大规模场景的经验**
在大规模场景中如果使用GUI会大量增加显存需求，可能导致程序崩溃。因此推荐先在小规模场景上试验，然后将几何尺寸等比例放大（保留粒子半径等模拟参数不变）。然后在命令行中使用"./bin/SPHSimulator.exe –-no-gui"选项运行程序，同时要在json中指定stopAt。