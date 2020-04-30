Outputter.cpp
629-650，新单元动力学单元应力输出
703-706，新单元信息
885-906，新单元静力学单元应力输出
以上内容合并后需取消注释，并将单元名修改为对应单元名称
另外亲测，vtk_Output->OutputVTK()必须在Tec_Output->OutputTecplot(1)之前，否则会读错