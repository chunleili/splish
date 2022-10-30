#pragma once
#include "SPlisHSPlasH/Common.h"

//点击时获取鼠标位置
struct Interactive
{
    Vector3r mouse_pos;

    //a singleton method to get the object
    static Interactive& get_inter()
    {
        static Interactive inter;
        return inter;
    }
    
    void operation()
    {
        printf("mouse pos in Inter:(%.3f,\t%.3f,\t%.3f)\n", mouse_pos[0],mouse_pos[1],mouse_pos[2]);
    }

    //把mouse_pos从外界传递给Interactive内部
    void get_mouse_pos(const Vector3r& rhs)
    {
        mouse_pos[0] = rhs[0];
        mouse_pos[1] = rhs[1];
        mouse_pos[2] = rhs[2];
    }

    //获取刚体的控制权，也就是位置的指针。
    void set_rb_pos(Vector3r* rb_pos)
    {
        (*rb_pos) = mouse_pos; // 获取并设定位置
        
        std::cout<< "rb_pos: "<< (*rb_pos)<<"\n";
    }
};