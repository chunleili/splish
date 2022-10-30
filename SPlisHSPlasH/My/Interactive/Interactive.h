#pragma once
#include "SPlisHSPlasH/Common.h"
// #include "GUI/OpenGL/MiniGL.h"
#include "Simulator/GUI/imgui/Simulator_GUI_imgui.h"
using namespace SPH;

//点击时获取鼠标位置
struct Interactive
{
    Vector3r mouse_pos;

    
    void operation()
    {
        printf("mouse pos in Inter:(%.3f,\t%.3f,\t%.3f)\n", mouse_pos[0],mouse_pos[1],mouse_pos[2]);
    }

    void set_mouse_pos(const Vector3r& rhs)
    {
        mouse_pos[0] = rhs[0];
        mouse_pos[1] = rhs[1];
        mouse_pos[2] = rhs[2];

    }
};