echo "                               "
echo "Configuring and building Thirdparty/DBoW2 ..."
cd Thirdparty/DBoW2
case $1 in 0) 
	mkdir build
	;;
esac
cd build
case $1 in 0) 
	cmake .. -DCMAKE_BUILD_TYPE=Release
	;;
esac
make -j




echo "                               "
echo "Configuring and building Thirdparty/g2o ..."
cd ../../g2o
case $1 in 0) 
	mkdir build
	;;
esac
cd build
case $1 in 0) 
	cmake .. -DCMAKE_BUILD_TYPE=Release
	;;
esac
make -j

cd ../../../




case $1 in 0) 
	echo "                               "
	echo "Uncompress vocabulary ..."
	cd Vocabulary
	tar -xf ORBvoc.txt.tar.gz
	cd ..
	;;
esac




echo "                               "
echo "Configuring and building ORB_SLAM2 ..."
case $1 in 0) 
	mkdir build
	;;
esac
cd build
case $1 in 0) 
	cmake .. -DCMAKE_BUILD_TYPE=Release
	;;
esac
make -j




echo "                               "
echo "Configuring and building ROS examples ..."
cd ../Examples/ROS/ORB_SLAM2
case $1 in 0) 
	mkdir build
	;;
esac
cd build
case $1 in 0) 
	cmake .. -DROS_BUILD_TYPE=Release
	;;
esac
make

echo "                               "

