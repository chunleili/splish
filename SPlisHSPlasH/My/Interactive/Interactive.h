#pragma once
#include "SPlisHSPlasH/Common.h"

//用户交互（键盘鼠标）的中介类，用于传递和处理数据
struct Interactive
{
    Vector3r mouse_pos;
    enum KEY{left=0, right, up, down, forward, backward};

    //a singleton method to get the object
    static Interactive& get_inter()
    {
        static Interactive inter;
        return inter;
    }
    
    //把mouse_pos从外界传递给Interactive内部
    void get_mouse_pos(const Vector3r& rhs)
    {
        mouse_pos[0] = rhs[0];
        mouse_pos[1] = rhs[1];
        mouse_pos[2] = rhs[2];
        printf("mouse pos in Inter:(%.3f,\t%.3f,\t%.3f)\n", mouse_pos[0],mouse_pos[1],mouse_pos[2]);
    }

    //获取键盘的输入：从GUI\OpenGL\MiniGL.cpp MiniGL::keyboard
    void get_key_input(enum KEY input)
    {
        if(input == KEY::down)
            printf("down!\n");
        else if(input == KEY::up)
            printf("up!\n");
        else if(input == KEY::left)
            printf("left!\n");
        else if(input == KEY::right)
            printf("righ!\n");
        else if(input == KEY::forward)
            printf("forward!\n");
        else if(input == KEY::backward)
            printf("backward!\n");
    }

    //获取刚体的控制权。
    void set_rb_pos(Vector3r& rb_pos)
    {
        (rb_pos) = mouse_pos; // 获取并设定位置
        
        std::cout<< "rb_pos: "<< (rb_pos)<<"\n";
        
        // rb_pos[2]=0;

        // if(rb_pos[2]<-0.5)
        //     rb_pos[2] = -.5;
        // if(rb_pos[2]>0.5)
        //     rb_pos[2] = .5;
    }
};