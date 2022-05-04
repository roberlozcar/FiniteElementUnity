-- Lua script.
p=tetview:new()
p:load_mesh("D:/Drive/TEACHING/Master/ASA FEM 2019/P4Enunciado/Assets/Resources/sphere_sim.1.ele")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
