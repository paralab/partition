#ifndef __DTYPES_H_
#define __DTYPES_H_

#include <mpi.h>
#include <complex>
#include <cstddef>

#include "../mesh-util/mesh-util.hpp"

/** 
  @file	dtypes.h
  @brief	Traits to determine MPI_DATATYPE from a C++ datatype
  @author	Hari Sundar, hsundar@gmail.com

  Traits to determine MPI_DATATYPE from a C++ datatype. For non standard
  C++ datatypes (like classes), we will need to define additional classes. An
  example is given for the case of the std. complex variable. Additional
  classes can be added as required.
 **/ 

namespace par {

/**
@class Mpi_datatype
@brief An abstract class used for communicating messages using user-defined datatypes.
 The user must implement the static member function "value()" that returns the MPI_Datatype corresponding 
to this user-defined datatype.
@author	Hari Sundar, hsundar@gmail.com
@see Mpi_datatype<bool>
*/
    template <typename T>
    class Mpi_datatype;

		template <typename T1, typename T2>
		class Mpi_pairtype;

#define HS_MPIDATATYPE(CTYPE, MPITYPE)	\
	template <>  \
		class Mpi_datatype<CTYPE> \
		{  \
		    public: \
			static MPI_Datatype value() {\
                          return MPITYPE;\
                        } \
		};

    HS_MPIDATATYPE(short,          MPI_SHORT)
    HS_MPIDATATYPE(int,            MPI_INT)
    HS_MPIDATATYPE(long,           MPI_LONG)
    HS_MPIDATATYPE(unsigned short, MPI_UNSIGNED_SHORT)
    HS_MPIDATATYPE(unsigned int,   MPI_UNSIGNED)
    HS_MPIDATATYPE(unsigned long,  MPI_UNSIGNED_LONG)
    HS_MPIDATATYPE(float,          MPI_FLOAT)
    HS_MPIDATATYPE(double,         MPI_DOUBLE)
    HS_MPIDATATYPE(long double,    MPI_LONG_DOUBLE)
    HS_MPIDATATYPE(long long,    MPI_LONG_LONG_INT)
    HS_MPIDATATYPE(char,    MPI_CHAR)
    HS_MPIDATATYPE(unsigned char,    MPI_UNSIGNED_CHAR)
    // HS_MPIDATATYPE(uint16_t,    MPI_UINT16_T)


    //PetscScalar is simply a typedef for double. Hence no need to explicitly
    //define an mpi_datatype for it.

#undef HS_MPIDATATYPE

#define HS_MPIPAIRDATATYPE(CTYPE1, CTYPE2, MPITYPE)	\
	template <>  \
		class Mpi_pairtype<CTYPE1, CTYPE2> \
		{  \
		    public: \
			static MPI_Datatype value() {\
                          return MPITYPE;\
                        } \
		};

    HS_MPIPAIRDATATYPE(float, 			int,		MPI_FLOAT_INT)
		HS_MPIPAIRDATATYPE(double, 			int, 		MPI_DOUBLE_INT)
		HS_MPIPAIRDATATYPE(long, 				int,    MPI_LONG_INT)		
	  HS_MPIPAIRDATATYPE(int, 				int,  	MPI_2INT)
		HS_MPIPAIRDATATYPE(short, 			int,  	MPI_SHORT_INT)
		HS_MPIPAIRDATATYPE(long double, int,    MPI_LONG_DOUBLE_INT)
				
#undef HS_MPIPAIRDATATYPE

		

    template <typename T>
    class Mpi_datatype<std::complex<T> > {
    public:
        static MPI_Datatype value()
        {
            static bool         first = true;
            static MPI_Datatype datatype;

            if (first) {
                first = false;
                MPI_Type_contiguous(2, Mpi_datatype<T>::value(), &datatype);
                MPI_Type_commit(&datatype);
            }

            return datatype;
        }
    };

    template <typename T1, typename T2>
    class Mpi_pairtype {
    public:
        static MPI_Datatype value()
        {
            static bool         first = true;
            static MPI_Datatype datatype;

            if (first) {
                first = false;
								int block[2];  
								MPI_Aint disp[2];
								MPI_Datatype type[2];
			
								block[0] = 1;
								block[0] = 1;
								type[0] = Mpi_datatype<T1>::value(); 
								type[1] = Mpi_datatype<T2>::value(); 
								disp[0] = 0; 
								disp[1] = sizeof(T1); 
								
								MPI_Type_create_struct(2, block, disp, type, &datatype);
                MPI_Type_commit(&datatype);
            }

            return datatype;
        }
    };

/**
@brief A template specialization of the abstract class Mpi_datatype. This can be used for communicating messages of type "bool"
@author Rahul Sampath, rahul.sampath@gmail.com
*/
    template <>
    class Mpi_datatype<bool> {
        static void bool_LOR(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for (int i = 0; i < (*len); i++) {
                ((static_cast<bool*>(inout))[i]) =
                  ( ((static_cast<bool*>(in))[i]) || 
                    ((static_cast<bool*>(inout))[i]) );
            }//end for	
        }//end function


        static void bool_LAND(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for (int i = 0; i < (*len); i++) {
                ((static_cast<bool*>(inout))[i]) =
                  ( ((static_cast<bool*>(in))[i]) && ((static_cast<bool*>(inout))[i]));
            }//end for	
        }//end function

    public: 
        /** 
          @brief User defined MPI_Operation to perform: second[i] = (first[i] && second[i]), 
          @remark first and second are 2 arrays of type bool and '&&' represents the 'Logical AND' operation. 
         **/
        static MPI_Op LAND() {
            static bool         first = true;
            static MPI_Op land;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<bool>::bool_LAND ,true ,&land);
            }

            return land;
        }

        /** 
          @brief User defined MPI_Operation to perform: second[i] = (first[i] || second[i]), 
          @remark first and second are 2 arrays of type bool and '||' represents the 'Logical OR' operation. 
         **/
        static MPI_Op LOR() {
            static bool         first = true;
            static MPI_Op lor;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<bool>::bool_LOR ,true ,&lor);
            }

            return lor;
        }

	  /** 
          @return the MPI_Datatype for the C++ datatype "bool"
         **/
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype datatype;
            if (first) {
                first = false;
                MPI_Type_contiguous(sizeof(bool), MPI_BYTE, &datatype);
                MPI_Type_commit(&datatype);
            }
            return datatype;
        }
    };

    template <>
    class Mpi_datatype<TetElementWithFaces> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "TetElementWithFaces"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[7] = {1, 1, 1, 1, 1, 1, 4};
                MPI_Datatype types[7] = {MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UINT64_T, MPI_UINT64_T};
                MPI_Aint offsets[7];
                offsets[0] = offsetof(TetElementWithFaces, element_tag);
                offsets[1] = offsetof(TetElementWithFaces, global_idx);
                offsets[2] = offsetof(TetElementWithFaces, x);
                offsets[3] = offsetof(TetElementWithFaces, y);
                offsets[4] = offsetof(TetElementWithFaces, z);
                offsets[5] = offsetof(TetElementWithFaces, morton_encoding);
                offsets[6] = offsetof(TetElementWithFaces, face_tags);


                MPI_Type_create_struct(7, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_datatype<HexElementWithFaces> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "HexElementWithFaces"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[7] = {1, 1, 1, 1, 1, 1, 6};
                MPI_Datatype types[7] = {MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UINT64_T, MPI_UINT64_T};
                MPI_Aint offsets[7];
                offsets[0] = offsetof(TetElementWithFaces, element_tag);
                offsets[1] = offsetof(TetElementWithFaces, global_idx);
                offsets[2] = offsetof(TetElementWithFaces, x);
                offsets[3] = offsetof(TetElementWithFaces, y);
                offsets[4] = offsetof(TetElementWithFaces, z);
                offsets[5] = offsetof(TetElementWithFaces, morton_encoding);
                offsets[6] = offsetof(TetElementWithFaces, face_tags);


                MPI_Type_create_struct(7, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_datatype<ElementWithFace> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "ElementWithFace"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[3] = {1, 1, 1};
                MPI_Datatype types[3] = {MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T};
                MPI_Aint offsets[3];
                offsets[0] = offsetof(ElementWithFace, element_tag);
                offsets[1] = offsetof(ElementWithFace, global_idx);
                offsets[2] = offsetof(ElementWithFace, face_tag);


                MPI_Type_create_struct(3, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_datatype<ElementWithTag> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "ElementWithTag"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[2] = {1, 1};
                MPI_Datatype types[2] = {MPI_UINT64_T, MPI_UINT64_T};
                MPI_Aint offsets[2];
                offsets[0] = offsetof(ElementWithTag, element_tag);
                offsets[1] = offsetof(ElementWithTag, global_idx);



                MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_datatype<ElementWithCoord> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "ElementWithCoord"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[5] = {1, 1, 1, 1, 1};
                MPI_Datatype types[5] = {MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
                MPI_Aint offsets[5];
                offsets[0] = offsetof(ElementWithCoord, element_tag);
                offsets[1] = offsetof(ElementWithCoord, global_idx);
                offsets[2] = offsetof(ElementWithCoord, x);
                offsets[3] = offsetof(ElementWithCoord, y);
                offsets[4] = offsetof(ElementWithCoord, z);




                MPI_Type_create_struct(5, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_pairtype<ElementWithTag,ElementWithTag> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "std::pair<ElementWithTag,ElementWithTag>"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                
                first = false;
                MPI_Datatype inner_type = Mpi_datatype<ElementWithTag>::value();
                
                int second_value_offset;
                MPI_Type_size(inner_type, &second_value_offset);
                int block_lengths[2] = {1, 1};
                MPI_Datatype types[2] = {inner_type, inner_type};
                MPI_Aint offsets[2];
                offsets[0] = 0;
                offsets[1] = static_cast<MPI_Aint>(second_value_offset);


                MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

} //end namespace

#endif

