# My
This is a forked repo from https://github.com/InteractiveComputerGraphics/SPlisHSPlasH

we add these models to the splishsplash:
1. Diffusion: see temperatureDiffusion
2. Local dynamic viscosity: see Viscosity/Viscosity_Casson
3. Plasticity: see My/Plasticity
4. ShapeMatching: see RigidBody
5. NonNewton: see NonNewton

这是对Jan Bender原本repo的fork。原repo为：https://github.com/InteractiveComputerGraphics/SPlisHSPlasH

所有的增改都尽量保证开闭原则，即宁增加，勿修改。尽量保证不动原有代码。凡是要增加新功能，都要新建一个类/自由函数，而不破坏原有代码。即使要修改原代码，也建议通过复制粘贴原代码，新建一个类，并在新的文件中更改。这样的好处是尽量让修改是可退的，保证出现差错之后，我们还能恢复到可以运行的版本。

我们对源码的修改都尽量几种放在SPlisHSPlasH/My这个文件夹下面。我们自己新增的场景都放在 data/MyScenes下面


我们在此基础上进行了许多增改：
1. Diffusion: 按照温度进行扩散。请看 temperatureDiffusion和Coagulation。
2. Viscosity_Casson: 请看 Viscosity/Viscosity_Casson。根据Casson公式计算黏度，并且根据温度动态调整黏度。
3. Plasticity: 弹塑性，请看塑性冲击的例子plastic_strike.json。请看 My/Plasticity。
4. (有BUG)ShapeMatching: 刚体。请看RigidBody_Bunny.json的例子。
5. (未完成)NonNewton: 请看 NonNewton。非牛顿剪切变稀和剪切变稠。
6. MyPartioReader：粒子的导入器，可以直接导入Houdini的 .bhclassic 文件(.bhclassic是Houdini旧版（12以前）的bgeo文件）。
7. 用户交互：现在可以用鼠标右键旋转视角，鼠标滚轮缩放视角，鼠标中键平移视角。更加符合人类直觉。
8. Interactive: 交互类。可以获取点击的鼠标世界位置（打印到屏幕上）。还可以控制刚体（通过control+WASDF控制第0号刚体）。
9. SurfaceParticles: 自由函数findSurfaceParticles。用于找到所有表面粒子（定义为邻居数量<10的粒子）
10.  Exporter_xyz: 导出粒子信息为xyz格式（就是每行输出ascii）。目前位于Simulator/exporter中
11.  Exporter_MyPartio: 导出的bgeo格式可以直接被Houdini读取。修复了原有partio导出器的BUG。目前位于Simulator/exporter中。
12. (未完成)MyTimeStep: 留空，以后将DFSPH搬到这里修改。
