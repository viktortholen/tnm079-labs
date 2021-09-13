# Copyright Link�ping University 2010
# Author Stefan Lindholm <stefan.lindholm@liu.se>

if(BUILD_LAB1)
	set(MY_MATH ${MY_MATH}
		Math/ConjugateGradient.h
		Math/ConstantVectorField.h
		Math/Function3D.h
		Math/TrilinearInterpolator.cpp
		Math/TrilinearInterpolator.h
		Math/Volume.h
		Math/VortexVectorField.h
	)
endif(BUILD_LAB1)
