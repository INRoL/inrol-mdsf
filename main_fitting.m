close all;
clear;
clc;

addpath("src/");
addpath("src/singlemesh");
addpath("src/optimizer/");
addpath("src/objparser/");

%% Object path
path = "./obj/fitting/fitting_test.obj";

%----- Output .obj path
output_path = "./res/fitting_res.obj";
% output_path = ""; % (No save)


%% Parameters initialization
para = initializer(20, 40);

%% Shape optimization
para.maxiter = 100;
para.maxfunc = 1000;
para.verbose = true;
para.showfig = true;
q = para.quat_id;

obj = mesh2dsf(path, q, para, output_path);

