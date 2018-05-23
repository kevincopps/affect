#!/usr/bin/env python

from libc.stddef cimport size_t
from libc.stdint cimport int64_t
from libc.stdint cimport int32_t
from libc.stdint cimport uint32_t


cdef extern from "omp.h":
    int omp_get_num_threads()
    void omp_set_num_threads(int)
    void omp_set_dynamic(int)

cdef extern from "array_util.h" nogil:
    void* ex_aligned_allocate(size_t size)
    void ex_aligned_copy_stride(double* source, size_t source_length, double* destination, int destination_stride)
    void ex_aligned_to_zero_based_int32(int32_t* array, size_t n)
    void ex_aligned_to_zero_based_int64(int64_t* array, size_t n)
    void ex_aligned_to_one_based_int32(int32_t* array, size_t n)
    void ex_aligned_to_one_based_int64(int64_t* array, size_t n)
    void ex_aligned_zero_uint32(uint32_t* a, size_t n)
    
cdef extern from "global_to_local.h" nogil:    
    uint32_t compute_global_to_local64(
        size_t num_entry,
        const int64_t * element_to_vertex_global,
        size_t * max_global_index,
        size_t * min_global_index,
        uint32_t * element_to_vertex_local,
        uint32_t * global_to_local)

    void fill_local_to_global64(
        size_t max_global_index,
        size_t min_global_index,
        const uint32_t * global_to_local,
        int64_t * local_to_global)

    uint32_t compute_global_to_local32(
        size_t num_entry,
        const int32_t * element_to_vertex_global,
        size_t * max_global_index,
        size_t * min_global_index,
        uint32_t * element_to_vertex_local,
        uint32_t * global_to_local)

    void fill_local_to_global32(
        size_t max_global_index,
        size_t min_global_index,
        const uint32_t * global_to_local,
        int32_t * local_to_global)
    
    
cdef extern from "exodusII.h":
    
    float API_VERS "EX_API_VERS"
    enum: EX_API_VERS_NODOT
    
    enum: EX_READ          # ex_open(): open file for reading(default)
    enum: EX_WRITE         # ex_open(): open existing file for appending.
    
    enum: EX_NOCLOBBER     # Don't overwrite existing database, default
    enum: EX_CLOBBER       # Overwrite existing database if it exists
    enum: EX_NORMAL_MODEL  # disable mods that permit storage of larger models
    enum: EX_LARGE_MODEL   # enable mods that permit storage of larger models
    enum: EX_NETCDF4       # use the hdf5-based netcdf4 output
    enum: EX_NOSHARE       # Do not open netcdf file in "share" mode
    enum: EX_SHARE         # Do open netcdf file in "share" mode
    enum: EX_NOCLASSIC     # Do not force netcdf to classic mode in netcdf4 mode
    
    enum:  EX_MAPS_INT64_DB
    enum:  EX_IDS_INT64_DB
    enum:  EX_BULK_INT64_DB
    enum:  EX_ALL_INT64_DB

    enum:  EX_MAPS_INT64_API
    enum:  EX_IDS_INT64_API
    enum:  EX_BULK_INT64_API
    enum:  EX_INQ_INT64_API
    enum:  EX_ALL_INT64_API

    enum:  EX_MPIIO
    enum:  EX_MPIPOSIX
    enum:  EX_PNETCDF
    
    ctypedef enum ex_inquiry:
        EX_INQ_FILE_TYPE   
        EX_INQ_API_VERS    
        EX_INQ_DB_VERS     
        EX_INQ_TITLE       
        EX_INQ_DIM         
        EX_INQ_NODES       
        EX_INQ_ELEM        
        EX_INQ_ELEM_BLK    
        EX_INQ_NODE_SETS   
        EX_INQ_NS_NODE_LEN 
        EX_INQ_SIDE_SETS   
        EX_INQ_SS_NODE_LEN 
        EX_INQ_SS_ELEM_LEN 
        EX_INQ_QA          
        EX_INQ_INFO        
        EX_INQ_TIME        
        EX_INQ_EB_PROP     
        EX_INQ_NS_PROP     
        EX_INQ_SS_PROP     
        EX_INQ_NS_DF_LEN   
        EX_INQ_SS_DF_LEN   
        EX_INQ_LIB_VERS    
        EX_INQ_EM_PROP     
        EX_INQ_NM_PROP     
        EX_INQ_ELEM_MAP    
        EX_INQ_NODE_MAP    
        EX_INQ_EDGE        
        EX_INQ_EDGE_BLK    
        EX_INQ_EDGE_SETS   
        EX_INQ_ES_LEN      
        EX_INQ_ES_DF_LEN   
        EX_INQ_EDGE_PROP   
        EX_INQ_ES_PROP     
        EX_INQ_FACE        
        EX_INQ_FACE_BLK    
        EX_INQ_FACE_SETS   
        EX_INQ_FS_LEN      
        EX_INQ_FS_DF_LEN   
        EX_INQ_FACE_PROP   
        EX_INQ_FS_PROP     
        EX_INQ_ELEM_SETS   
        EX_INQ_ELS_LEN     
        EX_INQ_ELS_DF_LEN  
        EX_INQ_ELS_PROP    
        EX_INQ_EDGE_MAP    
        EX_INQ_FACE_MAP    
        EX_INQ_COORD_FRAMES
        EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH
        EX_INQ_DB_MAX_USED_NAME_LENGTH
        EX_INQ_MAX_READ_NAME_LENGTH
        EX_INQ_DB_FLOAT_SIZE
        EX_INQ_NUM_CHILD_GROUPS
        EX_INQ_GROUP_PARENT
        EX_INQ_GROUP_ROOT
        EX_INQ_GROUP_NAME_LEN
        EX_INQ_GROUP_NAME
        EX_INQ_FULL_GROUP_NAME_LEN
        EX_INQ_FULL_GROUP_NAME
        EX_INQ_INVALID
    
    ctypedef enum ex_option_type:
        EX_OPT_MAX_NAME_LENGTH
        EX_OPT_COMPRESSION_TYPE
        EX_OPT_COMPRESSION_LEVEL
        EX_OPT_COMPRESSION_SHUFFLE
        EX_OPT_INTEGER_SIZE_API
        EX_OPT_INTEGER_SIZE_DB
  
  
    ctypedef enum ex_entity_type:
        EX_NODAL
        EX_NODE_BLOCK
        EX_NODE_SET
        EX_EDGE_BLOCK
        EX_EDGE_SET
        EX_FACE_BLOCK
        EX_FACE_SET
        EX_ELEM_BLOCK
        EX_ELEM_SET
        EX_SIDE_SET
        EX_ELEM_MAP
        EX_NODE_MAP
        EX_EDGE_MAP
        EX_FACE_MAP
        EX_GLOBAL
        EX_COORDINATE
        EX_INVALID
  
    ctypedef enum ex_options:
        EX_DEFAULT
        EX_VERBOSE
        EX_DEBUG
        EX_ABORT
        EX_NULLVERBOSE
  
    enum: EX_INVALID_ID

    enum: MAX_STR_LENGTH
    enum: MAX_NAME_LENGTH
    enum: MAX_LINE_LENGTH  # Maximum length of the database title or an information record
    enum: MAX_ERR_LENGTH
    
    int exerrval # shared error return value
    int exoptval # error reporting flag (default is quiet)

    # exerrval return values    
    enum: EX_MEMFAIL
    enum: EX_BADFILEMODE
    enum: EX_BADFILEID
    enum: EX_WRONGFILETYPE
    enum: EX_LOOKUPFAIL
    enum: EX_BADPARAM
    enum: EX_MSG
    enum: EX_PRTLASTMSG
    enum: EX_NOTROOTID
    enum: EX_NULLENTITY

    #
    # Type definitions
    #
    ctypedef int64_t ex_entity_id
  
    # The mechanism for passing double/float and int/int64_t both use a
    # void*; to avoid some confusion as to whether a function takes an
    # integer or a float/double, the following typedef is used for the
    # integer argument
    ctypedef void void_int
    
    #
    # Data structures
    #
    cdef struct ex_init_params:
        char title[MAX_LINE_LENGTH + 1]
        int64_t num_dim
        int64_t num_nodes
        int64_t num_edge
        int64_t num_edge_blk
        int64_t num_face
        int64_t num_face_blk
        int64_t num_elem
        int64_t num_elem_blk
        int64_t num_node_sets
        int64_t num_edge_sets
        int64_t num_face_sets
        int64_t num_side_sets
        int64_t num_elem_sets
        int64_t num_node_maps
        int64_t num_edge_maps
        int64_t num_face_maps
        int64_t num_elem_maps

    cdef struct ex_block:
        int64_t id
        ex_entity_type type
        char topology[MAX_STR_LENGTH+1]
        int64_t num_entry
        int64_t num_nodes_per_entry
        int64_t num_edges_per_entry
        int64_t num_faces_per_entry
        int64_t num_attribute

    cdef struct ex_set:
        int64_t id
        ex_entity_type type
        int64_t num_entry
        int64_t num_distribution_factor
        void_int* entry_list
        void_int* extra_list
        void* distribution_factor_list

    cdef struct ex_block_params:
        void_int*  edge_blk_id
        char** edge_type
        int* num_edge_this_blk
        int* num_nodes_per_edge
        int* num_attr_edge
        void_int* face_blk_id
        char** face_type
        int* num_face_this_blk
        int* num_nodes_per_face
        int* num_attr_face
        void_int* elem_blk_id
        char** elem_type
        int* num_elem_this_blk
        int* num_nodes_per_elem
        int* num_edges_per_elem
        int* num_faces_per_elem
        int* num_attr_elem
        int define_maps

    cdef struct ex_set_specs:
        void_int* sets_ids
        void_int* num_entries_per_set
        void_int* num_dist_per_set
        void_int* sets_entry_index
        void_int* sets_dist_index
        void_int* sets_entry_list
        void_int* sets_extra_list
        void* sets_dist_fact

    cdef struct ex_var_params:
        int num_glob
        int num_node
        int num_edge
        int num_face
        int num_elem
        int num_nset
        int num_eset
        int num_fset
        int num_sset
        int num_elset
        int* edge_var_tab
        int* face_var_tab
        int* elem_var_tab
        int* nset_var_tab
        int* eset_var_tab
        int* fset_var_tab
        int* sset_var_tab
        int* elset_var_tab

#
# Main Exodus file functions
#
   
    int ex_close(int exoid)

    int ex_copy(int in_exoid, int out_exoid)

    int ex_create(path, mode, comp_ws, io_ws)

    int ex_create_int(const char* path,
                      int cmode,
                      int* comp_ws,
                      int* io_ws,
                      int my_version)

    int ex_create_group(int parent_id,
                        const char* group_name)

    int ex_get_group_id(int exoid,
                        const char* group_name,
                        int* group_id)

    int ex_get_group_ids(int exoid,
                         int* num_children,
                         int* child_ids)

    int ex_get_all_times(int exoid,
                         void* time_values)

    int ex_get_coord_names(int exoid, char** coord_names)

    int ex_get_coord(int exoid,
                     void* x_coor,
                     void* y_coor,
                     void* z_coor)

    int ex_get_n_coord(int exoid,
                       int64_t start_node_num,
                       int64_t num_nodes,
                       void* x_coor,
                       void* y_coor,
                       void* z_coor)

    int ex_get_partial_coord(int exoid,
                             int64_t start_node_num,
                             int64_t num_nodes,
                             void* x_coor,
                             void* y_coor,
                             void* z_coor)

    int ex_get_ids(int exoid,
                   ex_entity_type obj_type,
                   void_int* ids)

    int ex_get_coordinate_frames(int exoid,
                                 int* nframes,
                                 void_int* cf_ids,
                                 void* pt_coordinates,
                                 char* tags)

    int ex_get_glob_vars(int exoid,
                         int time_step,
                         int num_glob_vars,
                         void* glob_var_vals)

    int ex_open_int(const char* path,
                    int mode,
                    int* comp_ws,
                    int* io_ws,
                    float* version,
                    int my_version)
    
    int ex_get_glob_var_time(int exoid,
                             int glob_var_index,
                             int beg_time_step,
                             int end_time_step,
                             void* glob_var_vals)

    int ex_get_info(int exoid, char** info)
    
    int ex_put_init_ext(int exoid,
                        const ex_init_params *param)
    
    int ex_get_init_ext(int exoid,
                        ex_init_params *param)
    
    int ex_get_init(int exoid,
                    char* title,
                    void_int* num_dim,
                    void_int* num_nodes,
                    void_int* num_elem,
                    void_int* num_elem_blk,
                    void_int* num_node_sets,
                    void_int* num_side_sets)

    int ex_put_init(int exoid,
                    const char* title,
                    int64_t num_dim,
                    int64_t num_nodes,
                    int64_t num_elem,
                    int64_t num_elem_blk,
                    int64_t num_node_sets,
                    int64_t num_side_sets)

    int ex_get_map_param(int exoid,
                         int* num_node_maps,
                         int* num_elem_maps)

    int ex_get_name(int exoid,
                    ex_entity_type obj_type,
                    ex_entity_id entity_id,
                    char* name)

    int ex_get_names(int exoid,
                     ex_entity_type obj_type,
                     char** names)

    int ex_get_nset_var_tab(int exoid,
                            int num_nodesets,
                            int num_nset_var,
                            int* nset_var_tab)

    int ex_get_n_nodal_var(int exoid,
                           int time_step,
                           int nodal_var_index,
                           int64_t start_node,
                           int64_t num_nodes,
                           void* nodal_var_vals)

    int ex_get_partial_nodal_var(int exoid,
                                 int time_step,
                                 int nodal_var_index,
                                 int64_t start_node,
                                 int64_t num_nodes,
                                 void* nodal_var_vals)

    int ex_get_prop_array(int exoid,
                          ex_entity_type obj_type,
                          const char* prop_name,
                          void_int* values)

    int ex_get_prop(int exoid,
                    ex_entity_type obj_type,
                    ex_entity_id obj_id,
                    const char* prop_name,
                    void_int* value)

    int ex_get_partial_num_map(int exoid,
                               ex_entity_type map_type,
                               ex_entity_id map_id,
                               int64_t ent_start,
                               int64_t ent_count,
                               void_int* elem_map)

    int ex_get_prop_names(int exoid,
                          ex_entity_type obj_type,
                          char** prop_names)

    int ex_get_qa(int exoid,
                  char* qa_record[][4])

    int ex_get_time(int exoid,
                    int time_step,
                    void* time_value)

    int ex_get_variable_names(int exoid,
                              ex_entity_type obj_type,
                              int num_vars,
                              char* var_names[])
    
    int ex_get_variable_name(int exoid,
                             ex_entity_type obj_type,
                             int var_num,
                             char* var_name)

    int ex_get_variable_param(int exoid,
                              ex_entity_type obj_type,
                              int* num_vars)

    int ex_get_object_truth_vector(int exoid,
                                   ex_entity_type var_type,
                                   ex_entity_id object_id,
                                   int num_var,
                                   int* var_vector)

    int ex_get_truth_table(int exoid,
                           ex_entity_type obj_type,
                           int num_blk,
                           int num_var,
                           int* var_tab)

    # ex_open() is implemented in exodus.pyx, as a call to ex_open_int
    #define ex_open(path, mode, comp_ws, io_ws, version) ex_open_int(path, mode, comp_ws, io_ws, version, EX_API_VERS_NODOT)

    int ex_open_int(const char* path,
                    int mode,
                    int* comp_ws,
                    int* io_ws,
                    float* version,
                    int my_version)

    int ex_put_attr_param(int exoid,
                          ex_entity_type obj_type,
                          ex_entity_id obj_id,
                          int num_attrs)

    int ex_get_attr_param(int exoid,
                          ex_entity_type obj_type,
                          ex_entity_id obj_id,
                          int* num_attrs)

    int ex_put_all_var_param(int exoid,
                             int num_g, int num_n,
                             int num_e, int* elem_var_tab,
                             int num_m, int* nset_var_tab,
                             int num_s, int* sset_var_tab)

    int ex_put_concat_elem_block(int exoid,
                                 const void_int* elem_blk_id,
                                 char* elem_type[],
                                 const void_int* num_elem_this_blk,
                                 const void_int* num_nodes_per_elem,
                                 const void_int* num_attr,
                                 int define_maps)

    int ex_put_coord_names(int exoid,
                           char* coord_names[])

    int ex_put_coord(int exoid,
                     const void* x_coor,
                     const void* y_coor,
                     const void* z_coor)

    int ex_put_n_coord(int exoid,
                       int64_t start_node_num,
                       int64_t num_nodes,
                       const void* x_coor,
                       const void* y_coor,
                       const void* z_coor)

    int ex_put_partial_coord(int exoid,
                             int64_t start_node_num,
                             int64_t num_nodes,
                             const void* x_coor,
                             const void* y_coor,
                             const void* z_coor)

    int ex_put_id_map(int exoid,
                      ex_entity_type obj_type,
                      const void_int* the_map)

    int ex_put_partial_id_map(int exoid,
                              ex_entity_type obj_type,
                              int64_t start_entity_num,
                              int64_t num_entities,
                              const void_int* the_map)

    int ex_put_n_elem_num_map(int exoid,
                              int64_t start_ent,
                              int64_t num_ents,
                              const void_int* the_map)

    int ex_put_n_node_num_map(int exoid,
                              int64_t start_ent,
                              int64_t num_ents,
                              const void_int* the_map)

    int ex_put_partial_elem_num_map(int exoid,
                                    int64_t start_ent,
                                    int64_t num_ents,
                                    const void_int* the_map)

    int ex_put_partial_node_num_map(int exoid,
                                    int64_t start_ent,
                                    int64_t num_ents,
                                    const void_int* the_map)

    int ex_get_id_map(int exoid,
                      ex_entity_type obj_type,
                      void_int* the_map)

    int ex_get_partial_id_map(int exoid,
                              ex_entity_type map_type,
                              int64_t start_entity_num,
                              int64_t num_entities,
                              void_int* id_map)

    int ex_put_coordinate_frames(int exoid,
                                 int nframes,
                                 const void_int* cf_ids,
                                 void* pt_coordinates,
                                 const char* tags)

    int ex_put_info(int exoid,
                    int num_info,
                    char* info[])

    int ex_put_map_param(int exoid,
                         int num_node_maps,
                         int num_elem_maps)

    int ex_put_name(int exoid,
                    ex_entity_type obj_type,
                    ex_entity_id entity_id,
                    const char* name)

    int ex_put_names(int exoid,
                     ex_entity_type obj_type,
                     char* names[])

    int ex_put_n_one_attr(int exoid,
                          ex_entity_type obj_type,
                          ex_entity_id obj_id,
                          int64_t start_num,
                          int64_t num_ent,
                          int attrib_index,
                          const void* attrib)

    int ex_put_partial_one_attr(int exoid,
                                ex_entity_type obj_type,
                                ex_entity_id obj_id,
                                int64_t start_num,
                                int64_t num_ent,
                                int attrib_index,
                                const void* attrib)

    int ex_put_prop(int exoid,
                    ex_entity_type obj_type,
                    ex_entity_id obj_id,
                    const char* prop_name,
                    ex_entity_id value)

    int ex_put_prop_array(int exoid,
                          ex_entity_type obj_type,
                          const char* prop_name,
                          const void_int* values)

    int ex_put_prop_names(int exoid,
                          ex_entity_type obj_type,
                          int num_props,
                          char** prop_names)

    int ex_put_qa(int exoid,
                  int num_qa_records,
                  char* qa_record[][4])

    int ex_put_time(int exoid,
                    int time_step,
                    const void* time_value)

    int ex_put_variable_name(int exoid,
                             ex_entity_type obj_type,
                             int var_num,
                             const char* var_name)

    int ex_put_variable_names(int exoid,
                              ex_entity_type obj_type,
                              int num_vars,
                              char* var_names[])

    int ex_put_variable_param(int exoid,
                              ex_entity_type obj_type,
                              int num_vars)

    int ex_put_truth_table(int exoid,
                           ex_entity_type obj_type,
                           int num_blk,
                           int num_var,
                           int* var_tab)

    int ex_update(int exoid)
    
    int ex_get_num_props(int exoid, ex_entity_type obj_type)
    
    int ex_large_model(int exoid)
    
    size_t ex_header_size(int exoid)

    void ex_err(const char* module_name, const char* message, int err_num)
    
    void ex_get_err(const char** msg, const char** func, int* errcode)
    
    int ex_opts(int options)
    
    int ex_inquire(int exoid, int inquiry, void_int*, float*, char*)
    
    int64_t ex_inquire_int(int exoid, int inquiry)
    
    int ex_int64_status(int exoid)
    
    int ex_set_int64_status(int exoid, int mode)

    # Note that the max name length setting is global at this time not specific
    # to a particular database however, the exoid option is passed to give
    # flexibility in the future to implement this on a database-by-database basis.
    int ex_set_max_name_length(int exoid, int length)

    int ex_set_option(int exoid, ex_option_type option, int option_value)

    # Write Node Edge Face or Element Number Map
    int ex_put_num_map(int exoid,
                       ex_entity_type map_type,
                       ex_entity_id map_id,
                       const void_int* the_map)

    #  Read Number Map
    int ex_get_num_map(int exoid,
                       ex_entity_type map_type,
                       ex_entity_id map_id,
                       void_int* the_map)

    #  Write Edge Face or Element Block Parameters
    int ex_put_block(int exoid,
                     ex_entity_type blk_type,
                     ex_entity_id blk_id,
                     const char* entry_descrip,
                     int64_t num_entries_this_blk,
                     int64_t num_nodes_per_entry,
                     int64_t num_edges_per_entry,
                     int64_t num_faces_per_entry,
                     int64_t num_attr_per_entry)

    # deprecated @see ex_get_block_param()
    int ex_get_block(int exoid,
                     ex_entity_type blk_type,
                     ex_entity_id blk_id,
                     char* elem_type,
                     void_int* num_entries_this_blk,
                     void_int* num_nodes_per_entry,
                     void_int* num_edges_per_entry,
                     void_int* num_faces_per_entry,
                     void_int* num_attr_per_entry)

    #  Read Edge Face or Element Block Parameters
    int ex_get_block_param(int exoid,
                           ex_block *block)

    int ex_put_block_param(int exoid,
                           const ex_block block)

    #  Write All Edge Face and Element Block Parameters
    int ex_put_concat_all_blocks(int exoid,
                                 const ex_block_params *param)

    int ex_put_entity_count_per_polyhedra(int exoid,
                                          ex_entity_type blk_type,
                                          ex_entity_id blk_id,
                                          const int* entity_counts)

    int ex_get_entity_count_per_polyhedra(int exoid,
                                          ex_entity_type blk_type,
                                          ex_entity_id blk_id,
                                          int* entity_counts)

    #  Write Edge Face or Element Block Connectivity
    int ex_put_conn(int exoid,
                    ex_entity_type blk_type,
                    ex_entity_id blk_id,
                    const void_int* node_conn,
                    const void_int* elem_edge_conn,
                    const void_int* elem_face_conn)

    #  Read Edge Face or Element Block Connectivity
    int ex_get_conn(int exoid,
                    ex_entity_type blk_type,
                    ex_entity_id blk_id,
                    void_int* nodeconn,
                    void_int* edgeconn,
                    void_int* faceconn)

    #  Read Partial Edge Face or Element Block Connectivity
    int ex_get_n_conn(int exoid,
                      ex_entity_type blk_type,
                      ex_entity_id blk_id,
                      int64_t start_num,
                      int64_t num_ent,
                      void_int* nodeconn,
                      void_int* edgeconn,
                      void_int* faceconn)

    int ex_get_partial_conn(int exoid,
                            ex_entity_type blk_type,
                            ex_entity_id blk_id,
                            int64_t start_num,
                            int64_t num_ent,
                            void_int* nodeconn,
                            void_int* edgeconn,
                            void_int* faceconn)

    #  Write Edge Face or Element Block Attributes
    int ex_put_attr(int exoid,
                      ex_entity_type blk_type,
                      ex_entity_id blk_id,
                      const void* attrib)

    int ex_put_partial_attr(int exoid,
                            ex_entity_type blk_type,
                            ex_entity_id blk_id,
                            int64_t start_entity,
                            int64_t num_entity,
                            const void* attrib)

    #  Read Edge Face or Element Block Attributes
    int ex_get_attr(int exoid,
                    ex_entity_type obj_type,
                    ex_entity_id obj_id,
                    void* attrib)

    int ex_get_n_attr(int exoid,
                      ex_entity_type obj_type,
                      ex_entity_id obj_id,
                      int64_t start_num,
                      int64_t num_ent,
                      void* attrib)

    int ex_get_partial_attr(int exoid,
                            ex_entity_type obj_type,
                            ex_entity_id obj_id,
                            int64_t start_num,
                            int64_t num_ent,
                            void* attrib)

    #  Write One Edge Face or Element Block Attribute
    int ex_put_one_attr(int exoid,
                        ex_entity_type obj_type,
                        ex_entity_id obj_id,
                        int attrib_index,
                        const void* attrib)

    #  Read One Edge Face or Element Block Attribute
    int ex_get_one_attr(int exoid,
                        ex_entity_type obj_type,
                        ex_entity_id obj_id,
                        int attrib_index,
                        void* attrib)

    #  Read One Edge Face or Element Block Attribute
    int ex_get_n_one_attr(int exoid,
                          ex_entity_type obj_type,
                          ex_entity_id obj_id,
                          int64_t start_num,
                          int64_t num_ent,
                          int attrib_index,
                          void* attrib)

    int ex_get_partial_one_attr(int exoid,
                                ex_entity_type obj_type,
                                ex_entity_id obj_id,
                                int64_t start_num,
                                int64_t num_ent,
                                int attrib_index,
                                void* attrib)

    #  Write Edge Face or Element Block Attribute Names
    int ex_put_attr_names(int exoid,
                          ex_entity_type blk_type,
                          ex_entity_id blk_id,
                          char** names)

    #  Read Edge Face or Element Block Attribute Names
    int ex_get_attr_names(int exoid,
                          ex_entity_type obj_type,
                          ex_entity_id obj_id,
                          char** names)

    #  Write Node Edge Face or Side Set Parameters
    int ex_put_set_param(int exoid,
                         ex_entity_type set_type,
                         ex_entity_id set_id,
                         int64_t num_entries_in_set,
                         int64_t num_dist_fact_in_set)

    #  Read Node Edge Face or Side Set Parameters
    int ex_get_set_param(int exoid,
                         ex_entity_type set_type,
                         ex_entity_id set_id,
                         void_int* num_entry_in_set,
                         void_int* num_dist_fact_in_set)

    #  Write a Node Edge Face or Side Set
    int ex_put_set(int exoid,
                   ex_entity_type set_type,
                   ex_entity_id set_id,
                   const void_int* set_entry_list,
                   const void_int* set_extra_list)

    int ex_put_partial_set(int exoid,
                           ex_entity_type set_type,
                           ex_entity_id set_id,
                           int64_t offset,
                           int64_t count,
                           const void_int* set_entry_list,
                           const void_int* set_extra_list)

    #  Read a Node Edge Face or Side Set
    int ex_get_set(int exoid,
                   ex_entity_type set_type,
                   ex_entity_id set_id,
                   void_int* set_entry_list,
                   void_int* set_extra_list)

    #  Write Node Edge Face or Side Set Distribution Factors
    int ex_put_set_dist_fact(int exoid,
                             ex_entity_type set_type,
                             ex_entity_id set_id,
                             const void* set_dist_fact)

    #  Read Node Edge Face or Side Set Distribution Factors
    int ex_get_set_dist_fact(int exoid,
                             ex_entity_type set_type,
                             ex_entity_id set_id,
                             void* set_dist_fact)

    int ex_get_partial_set_dist_fact(int exoid,
                                     ex_entity_type set_type,
                                     ex_entity_id set_id,
                                     int64_t offset,
                                     int64_t num_to_put,
                                     void* set_dist_fact)

    #  Write Concatenated Node Edge Face or Side Sets
    int ex_put_concat_sets(int exoid,
                           ex_entity_type set_type,
                           const ex_set_specs* set_specs)

    #  Read Concatenated Node Edge Face or Side Sets
    int ex_get_concat_sets(int exoid,
                           ex_entity_type set_type,
                           const ex_set_specs* set_specs)

    #  Write Concatenated Node Edge Face or Side Sets
    int ex_put_sets(int exoid,
                    size_t set_count,
                    const ex_set* sets)

    #  Read Concatenated Node Edge Face or Side Sets
    int ex_get_sets(int exoid,
                    size_t set_count,
                    ex_set* sets)

    # (MODIFIED) Write All Results Variables Parameters
    int ex_put_all_var_param_ext(int exoid,
                                 const ex_var_params* vp)

    #  Write Edge Face or Element Variable Values on Blocks or Sets at a Time Step
    int ex_put_var(int exoid,
                   int time_step,
                   ex_entity_type var_type,
                   int var_index,
                   ex_entity_id obj_id,
                   int64_t num_entries_this_obj,
                   const void* var_vals)

    #  Write Partial Edge Face or Element Variable Values on Blocks or Sets at a Time Step
    int ex_put_n_var(int exoid,
                     int time_step,
                     ex_entity_type var_type,
                     int var_index,
                     ex_entity_id obj_id,
                     int64_t start_index,
                     int64_t num_entities,
                     const void* var_vals)

    int ex_put_partial_var(int exoid,
                           int time_step,
                           ex_entity_type var_type,
                           int var_index,
                           ex_entity_id obj_id,
                           int64_t start_index,
                           int64_t num_entities,
                           const void* var_vals)

    #  Read Edge Face or Element Variable Values Defined On Blocks or Sets at a Time Step
    int ex_get_var(int exoid,
                   int time_step,
                   ex_entity_type var_type,
                   int var_index,
                   ex_entity_id obj_id,
                   int64_t num_entry_this_obj,
                   void* var_vals)

    #  Read Partial Edge Face or Element Variable Values on Blocks or Sets at a Time Step
    int ex_get_n_var(int exoid,
                     int time_step,
                     ex_entity_type var_type,
                     int var_index,
                     ex_entity_id obj_id,
                     int64_t start_index,
                     int64_t num_entities,
                     void* var_vals)

    int ex_get_n_elem_var(int exoid,
                          int time_step,
                          int elem_var_index,
                          ex_entity_id elem_blk_id,
                          int64_t num_elem_this_blk,
                          int64_t start_elem_num,
                          int64_t num_elem,
                          void* elem_var_vals)

    int ex_get_partial_var(int exoid,
                           int time_step,
                           ex_entity_type var_type,
                           int var_index,
                           ex_entity_id obj_id,
                           int64_t start_index,
                           int64_t num_entities,
                           void* var_vals)

    int ex_get_partial_elem_var(int exoid,
                                int time_step,
                                int elem_var_index,
                                ex_entity_id elem_blk_id,
                                int64_t num_elem_this_blk,
                                int64_t start_elem_num,
                                int64_t num_elem,
                                void* elem_var_vals)

    #  Read Edge Face or Element Variable Values Defined On Blocks or Sets Through Time
    int ex_get_var_time(int exoid,
                        ex_entity_type var_type,
                        int var_index,
                        int64_t varid,
                        int beg_time_step,
                        int end_time_step,
                        void* var_vals)

    int ex_cvt_nodes_to_sides(int exoid,
                              void_int* num_elem_per_set,
                              void_int* num_nodes_per_set,
                              void_int* side_sets_elem_index,
                              void_int* side_sets_node_index,
                              void_int* side_sets_elem_list,
                              void_int* side_sets_node_list,
                              void_int* side_sets_side_list)

    # Can be replaced by ex_put_var ...
    int ex_put_nodal_var(int exoid,
                         int time_step,
                         int nodal_var_index,
                         int64_t num_nodes,
                         const void* nodal_var_vals)

    int ex_put_n_nodal_var(int exoid,
                           int time_step,
                           int nodal_var_index,
                           int64_t start_node,
                           int64_t num_nodes,
                           const void* nodal_var_vals)

    int ex_put_partial_nodal_var(int exoid,
                                 int time_step,
                                 int nodal_var_index,
                                 int64_t start_node,
                                 int64_t num_nodes,
                                 const void* nodal_var_vals)

    int ex_get_partial_elem_map(int exoid,
                                ex_entity_id map_id,
                                int64_t ent_start,
                                int64_t ent_count,
                                void_int* elem_map)

    int ex_put_partial_elem_map(int exoid,
                                ex_entity_id map_id,
                                int64_t ent_start,
                                int64_t ent_count,
                                const void_int* elem_map)

    int ex_put_partial_num_map(int exoid,
                               ex_entity_type map_type,
                               ex_entity_id map_id,
                               int64_t ent_start,
                               int64_t ent_count,
                               const void_int* the_map)

    int ex_put_partial_set_dist_fact(int exoid,
                                     ex_entity_type set_type,
                                     ex_entity_id set_id,
                                     int64_t offset,
                                     int64_t num_to_put,
                                     const void* set_dist_fact)

    int ex_get_concat_side_set_node_count(int exoid,
                                          int* side_set_node_cnt_list)

    int ex_get_side_set_node_list_len(int exoid,
                                      ex_entity_id side_set_id,
                                      void_int* side_set_node_list_len)

    int ex_get_side_set_node_count(int exoid,
                                   ex_entity_id side_set_id,
                                   int* side_set_node_cnt_list)

    int ex_get_side_set_node_list(int exoid,
                                  ex_entity_id side_set_id,
                                  void_int* side_set_node_cnt_list,
                                  void_int* side_set_node_list)

    #---
    # Functions pulled from nemesis library and incorporated into exodus...
    #---
    
    #---
    # Initial Information Routines
    #---
    int ex_get_init_info(int exoid,             # NemesisI file ID
                         int* num_proc,         # Number of processors
                         int* num_proc_in_f,    # Number of procs in this file
                         char* ftype)

    int ex_put_init_info(int exoid,             # NemesisI file ID
                         int num_proc,          # Number of processors
                         int num_proc_in_f,     # Number of procs in this file
                         char* ftype)

    int ex_get_init_global(int exoid,                   # NemesisI file ID
                           void_int* num_nodes_g,       # Number of global FEM nodes
                           void_int* num_elems_g,       # Number of global FEM elements
                           void_int* num_elem_blks_g,   # Number of global elem blocks
                           void_int* num_node_sets_g,   # Number of global node sets
                           void_int* num_side_sets_g)   # Number of global side sets
    
    int ex_put_init_global(int exoid,                   # NemesisI file ID
                           int64_t num_nodes_g,	        # Number of global FEM nodes
                           int64_t num_elems_g,	        # Number of global FEM elements
                           int64_t num_elem_blks_g,	# Number of global elem blocks
                           int64_t num_node_sets_g,	# Number of global node sets
                           int64_t num_side_sets_g)	# Number of global side sets

    #---
    # Loadbalance Parameter Routines
    #---
    int ex_get_loadbal_param(int exoid,                 # NetCDF/Exodus file ID
                             void_int* num_int_nodes,   # Number of internal FEM nodes
                             void_int* num_bor_nodes,   # Number of border FEM nodes
                             void_int* num_ext_nodes,   # Number of external FEM nodes
                             void_int* num_int_elems,   # Number of internal FEM elems
                             void_int* num_bor_elems,   # Number of border FEM elems
                             void_int* num_node_cmaps,  # Number of nodal comm maps
                             void_int* num_elem_cmaps,  # Number of elemental comm maps
                             int processor)             # Processor ID

    int ex_put_loadbal_param(int exoid,                 # NemesisI file ID
                             int64_t num_int_nodes,     # Number of internal FEM nodes
                             int64_t num_bor_nodes,     # Number of border FEM nodes
                             int64_t num_ext_nodes,     # Number of external FEM nodes
                             int64_t num_int_elems,     # Number of internal FEM elems
                             int64_t num_bor_elems,     # Number of border FEM elems
                             int64_t num_node_cmaps,    # Number of nodal comm maps
                             int64_t num_elem_cmaps,    # Number of elemental comm maps
                             int processor)             # Processor ID


    int ex_put_loadbal_param_cc(int exoid,                  # NetCDF/Exodus file ID
                                void_int* num_int_nodes,    # Number of internal node IDs
                                void_int* num_bor_nodes,    # Number of border node IDs
                                void_int* num_ext_nodes,    # Number of external node IDs
                                void_int* num_int_elems,    # Number of internal elem IDs
                                void_int* num_bor_elems,    # Number of border elem IDs
                                void_int* num_node_cmaps,   # Number of nodal comm maps
                                void_int* num_elem_cmaps)   # Number of elem comm maps

    #---
    # NS, SS & EB Global Parameter Routines
    #---
    int ex_get_ns_param_global(int exoid,                   # NetCDF/Exodus file ID
                               void_int* ns_ids_glob,       # Global IDs of node sets
                               void_int* ns_n_cnt_glob,     # Count of nodes in node sets
                               void_int* ns_df_cnt_glob)    # Count of dist. factors in ns

    int ex_put_ns_param_global(int exoid,                   # NemesisI file ID
                               void_int* global_ids,	    # Vector of global node-set IDs
                               void_int* global_n_cnts,     # Vector of node counts in node-sets
                               void_int* global_df_cnts)    # Vector of dist factor counts in node-sets

    int ex_get_ss_param_global(int exoid,           # NetCDF/Exodus file ID
                       void_int* ss_ids_glob,       # Global side-set IDs
                       void_int* ss_s_cnt_glob,     # Global side count
                       void_int* ss_df_cnt_glob)    # Global dist. factor count

    int ex_put_ss_param_global(int exoid,           # NemesisI file ID
                       void_int* global_ids,        # Vector of global side-set IDs
                       void_int* global_el_cnts,    # Vector of element/side
                                                    # counts in each side set
                       void_int* global_df_cnts)    # Vector of dist. factor
                                                    # counts in each side set

    int ex_get_eb_info_global(int exoid,                # NemesisI file ID
                              void_int* el_blk_ids,     # Vector of global element IDs
                              void_int* el_blk_cnts)    # Vector of global element counts

    int ex_put_eb_info_global(int exoid,		# NemesisI file ID
                              void_int* el_blk_ids,	# Vector of global element IDs
                              void_int* el_blk_cnts)	# Vector of global element counts

    #---
    # NS, SS & EB Subset Routines
    #---
    int ex_get_n_side_set(int exoid,                            # NetCDF/Exodus file ID
                          ex_entity_id side_set_id,             # Side-set ID to read
                          int64_t start_side_num,               # Starting element number
                          int64_t num_sides,                    # Number of sides to read
                          void_int* side_set_elem_list,         # List of element IDs
                          void_int* side_set_side_list)         # List of side IDs

    int ex_put_n_side_set(int exoid,                            # NetCDF/Exodus file ID
                          ex_entity_id side_set_id,             # Side-set ID to write
                          int64_t start_side_num,               # Starting element number
                          int64_t num_sides,                    # Number of sides to write
                          const void_int* side_set_elem_list,   # List of element IDs
                          const void_int* side_set_side_list)   # List of side IDs

    int ex_get_n_side_set_df(int exoid,	                        # NetCDF/Exodus file ID
                             ex_entity_id side_set_id,          # Side-set ID
                             int64_t start_num,                 # Starting df number
                             int64_t num_df_to_get,             # Number of df's to read
                             void* side_set_df)                 # Distribution factors

    int ex_put_n_side_set_df(int exoid,                         # NetCDF/Exodus file ID
                             ex_entity_id side_set_id,          # Side-set ID
                             int64_t start_num,                 # Starting df number
                             int64_t num_df_to_get,             # Number of df's to write
                             void* side_set_df)                 # Distribution factors

    int ex_get_n_node_set(int exoid,                            # NetCDF/Exodus file ID
                          ex_entity_id  node_set_id,            # Node set ID
                          int64_t start_node_num,               # Node index to start reading at
                          int64_t num_node,                     # Number of nodes to read
                          void_int* node_set_node_list)         # List of nodes in node set

    int ex_put_n_node_set(int exoid,                            # NetCDF/Exodus file ID
                          ex_entity_id  node_set_id,            # Node set ID
                          int64_t start_node_num,               # Node index to start writing at
                          int64_t num_node,                     # Number of nodes to write
                          const void_int* node_set_node_list)   # List of nodes in node set

    int ex_get_n_node_set_df(int exoid,                         # NetCDF/Exodus file ID
                             ex_entity_id node_set_id,	        # Node-set ID
                             int64_t start_num,                 # Starting df number
                             int64_t num_df_to_get,             # Number of df's to read
                             void* node_set_df)                 # Distribution factors

    int ex_put_n_node_set_df(int exoid,                         # NetCDF/Exodus file ID
                             ex_entity_id node_set_id,          # Node-set ID
                             int64_t start_num,                 # Starting df number
                             int64_t num_df_to_get,             # Number of df's to write
                             void* node_set_df)                 # Distribution factors

    int ex_get_n_elem_conn(int exoid,                           # NetCDF/Exodus file ID
                           ex_entity_id elem_blk_id,            # Element block ID
                           int64_t start_elem_num,              # Starting position to read from
                           int64_t num_elems,	                # Number of elements to read
                           void_int* connect)                   # Connectivity vector

    int ex_put_n_elem_conn(int exoid,                           # NetCDF/Exodus file ID
                           ex_entity_id elem_blk_id,            # Element block ID
                           int64_t start_elem_num,              # Starting position to write to
                           int64_t num_elems,                   # Number of elements to write
                           const void_int* connect)             # Connectivity vector

    int ex_get_n_elem_attr(int exoid,                           # NetCDF/Exodus file ID
                           ex_entity_id elem_blk_id,            # Element block ID
                           int64_t start_elem_num,              # Starting position to read from
                           int64_t num_elems,                   # Number of elements to read
                           void* attrib)                        # Attribute

    int ex_put_n_elem_attr(int exoid,                           # NetCDF/Exodus file ID
                           ex_entity_id elem_blk_id,            # Element block ID
                           int64_t start_elem_num,              # Starting position to write to
                           int64_t num_elems,                   # Number of elements to write
                           void* attrib)                        # Attribute

    int ex_get_partial_side_set(int exoid,                      # NetCDF/Exodus file ID
                                ex_entity_id side_set_id,       # Side-set ID to read 
                                int64_t start_side_num,         # Starting element number
                                int64_t num_sides,              # Number of sides to read
                                void_int* side_set_elem_list,   # List of element IDs
                                void_int* side_set_side_list)   # List of side IDs

    int ex_put_partial_side_set(int exoid,                          # NetCDF/Exodus file ID
                                ex_entity_id side_set_id,           # Side-set ID to write
                                int64_t start_side_num,             # Starting element number
                                int64_t num_sides,                  # Number of sides to write
                                const void_int* side_set_elem_list, # List of element IDs
                                const void_int* side_set_side_list) # List of side IDs

    int ex_get_partial_side_set_df(int exoid,                   # NetCDF/Exodus file ID
                                   ex_entity_id side_set_id,    # Side-set ID
                                   int64_t start_num,           # Starting df number
                                   int64_t num_df_to_get,       # Number of df's to read
                                   void* side_set_df)           # Distribution factors

    int ex_put_partial_side_set_df(int exoid,                   # NetCDF/Exodus file ID
                                   ex_entity_id side_set_id,    # Side-set ID
                                   int64_t start_num,           # Starting df number
                                   int64_t num_df_to_get,       # Number of df's to write
                                   void* side_set_df)           # Distribution factors

    int ex_get_partial_node_set(int exoid,                      # NetCDF/Exodus file ID
                                ex_entity_id  node_set_id,      # Node set ID
                                int64_t start_node_num,         # Node index to start reading at
                                int64_t num_node,               # Number of nodes to read
                                void_int* node_set_node_list)   # List of nodes in node set

    int ex_put_partial_node_set(int exoid,	                        # NetCDF/Exodus file ID
                                ex_entity_id  node_set_id,          # Node set ID
                                int64_t start_node_num,             # Node index to start writing at
                                int64_t num_node,                   # Number of nodes to write
                                const void_int* node_set_node_list) # List of nodes in node set

    int ex_get_partial_node_set_df(int exoid,                   # NetCDF/Exodus file ID
                                   ex_entity_id node_set_id,    # Node-set ID
                                   int64_t start_num,           # Starting df number
                                   int64_t num_df_to_get,       # Number of df's to read
                                   void* node_set_df)           # Distribution factors

    int ex_put_partial_node_set_df(int exoid,                   # NetCDF/Exodus file ID
                                   ex_entity_id node_set_id,    # Node-set ID
                                   int64_t start_num,	        # Starting df number
                                   int64_t num_df_to_get,	    # Number of df's to write
                                   void* node_set_df)		    # Distribution factors

    int ex_get_partial_elem_conn(int exoid,                     # NetCDF/Exodus file ID
                                 ex_entity_id elem_blk_id,      # Element block ID
                                 int64_t start_elem_num,        # Starting position to read from
                                 int64_t num_elems,	            # Number of elements to read
                                 void_int* connect)	            # Connectivity vector

    int ex_put_partial_elem_conn(int exoid,                     # NetCDF/Exodus file ID
                                 ex_entity_id elem_blk_id,      # Element block ID
                                 int64_t start_elem_num,        # Starting position to write to
                                 int64_t num_elems,             # Number of elements to write
                                 const void_int* connect)       # Connectivity vector

    int ex_get_partial_elem_attr(int exoid,                     # NetCDF/Exodus file ID
                                 ex_entity_id elem_blk_id,      # Element block ID
                                 int64_t start_elem_num,        # Starting position to read from
                                 int64_t num_elems,             # Number of elements to read
                                 void* attrib)                  # Attribute

    int ex_put_partial_elem_attr(int exoid,                     # NetCDF/Exodus file ID
                                 ex_entity_id elem_blk_id,      # Element block ID
                                 int64_t start_elem_num,        # Starting position to write to
                                 int64_t num_elems,             # Number of elements to write
                                 void* attrib)                  # Attribute

    int ex_get_elem_type(int exoid,                             # NetCDF/Exodus file ID
                         ex_entity_id elem_blk_id,              # Element block ID
                         char* elem_type)                       # The name of the element type

    #---
    # Variable Routines
    #---

    int ex_put_elem_var_slab(int exoid,                 # NetCDF/Exodus file ID
                             int time_step,             # time index
                             int elem_var_index,        # elemental variable index
                             ex_entity_id elem_blk_id,  # elemental block id
                             int64_t start_pos,         # Starting position to write to
                             int64_t num_vals,          # Number of elements to write
                             void* elem_var_vals)       # variable values

    int ex_put_nodal_var_slab(int exoid,                # NetCDF/Exodus file ID
                              int time_step,            # The time step index
                              int nodal_var_index,      # Nodal variable index
                              int64_t start_pos,        # Start position for write
                              int64_t num_vals,         # Number of nodal variables
                              void* nodal_var_vals)     # Nodal variable values

    #---
    # Number Map Routines
    #---
    int ex_get_n_elem_num_map(int exoid,                # NetCDF/Exodus file ID
                              int64_t start_ent,        # Starting position to read from
                              int64_t num_ents,         # Number of elements to read
                              void_int* elem_map)       # element map numbers

    int ex_get_n_node_num_map(int exoid,                # NetCDF/Exodus file ID
                              int64_t start_ent,        # starting node number
                              int64_t num_ents,         # number of nodes to read
                              void_int* node_map)       # vector for node map

    int ex_get_partial_elem_num_map(int exoid,          # NetCDF/Exodus file ID
                                    int64_t start_ent,  # Starting position to read from
                                    int64_t num_ents,   # Number of elements to read
                                    void_int* elem_map) # element map numbers

    int ex_get_partial_node_num_map(int exoid,          # NetCDF/Exodus file ID
                                    int64_t start_ent,  # starting node number
                                    int64_t num_ents,   # number of nodes to read
                                    void_int* node_map) # vector for node map

    int ex_get_processor_node_maps(int exoid,           # NetCDF/Exodus file ID
                                   void_int* node_mapi,	# Internal FEM node IDs
                                   void_int* node_mapb,	# Border FEM node IDs
                                   void_int* node_mape,	# External FEM node IDs
                                   int processor)       # Processor IDs

    int ex_put_processor_node_maps(int exoid,           # NetCDF/Exodus file ID
                                   void_int* node_mapi,	# Internal FEM node IDs
                                   void_int* node_mapb,	# Border FEM node IDs
                                   void_int* node_mape,	# External FEM node IDs
                                   int processor)       # This processor ID

    int ex_get_processor_elem_maps(int exoid,		    # NetCDF/Exodus file ID
                                   void_int* elem_mapi,	# Internal element IDs
                                   void_int* elem_mapb,	# Border element IDs
                                   int processor)       # Processor ID

    int ex_put_processor_elem_maps(int exoid,		    # NetCDF/Exodus file ID
                                   void_int* elem_mapi,	# Internal FEM element IDs
                                   void_int* elem_mapb,	# Border FEM element IDs
                                   int processor)       # This processor ID

    #---
    # Communications Maps Routines
    #---

    int ex_get_cmap_params(int exoid,                       # NetCDF/Exodus file ID
                           void_int* node_cmap_ids,         # Nodal comm. map IDs
                           void_int* node_cmap_node_cnts,   # Number of nodes in each map
                           void_int* elem_cmap_ids,         # Elemental comm. map IDs
                           void_int* elem_cmap_elem_cnts,   # Number of elems in each map
                           int processor)                   # This processor ID

    int ex_put_cmap_params(int exoid,                       # NetCDF/Exodus file ID
                           void_int* node_map_ids,          # Node map IDs
                           void_int* node_map_node_cnts,    # Nodes in nodal comm
                           void_int* elem_map_ids,          # Elem map IDs
                           void_int* elem_map_elem_cnts,    # Elems in elemental comm
                           int64_t processor)               # This processor ID

    int ex_put_cmap_params_cc(int exoid,                    # NetCDF/Exodus file ID
                              void_int* node_map_ids,	    # Node map IDs
                              void_int* node_map_node_cnts, # Nodes in nodal comm
                              void_int* node_proc_ptrs,     # Pointer into array for
                                                            # node maps
                              void_int* elem_map_ids,	    # Elem map IDs
                              void_int* elem_map_elem_cnts, # Elems in elemental comm
                              void_int* elem_proc_ptrs)	    # Pointer into array for
                                                            # elem maps

    int ex_get_node_cmap(int exoid,             # NetCDF/Exodus file ID
                         ex_entity_id  map_id,  # Map ID
                         void_int* node_ids,    # FEM node IDs
                         void_int* proc_ids,    # Processor IDs
                         int processor)         # This processor ID

    int ex_put_node_cmap(int exoid,             # NetCDF/Exodus file ID
                         ex_entity_id  map_id,	# Nodal comm map ID
                         void_int* node_ids,	# FEM node IDs
                         void_int* proc_ids,    # Processor IDs
                         int processor)	        # This processor ID

    int ex_get_elem_cmap(int exoid,             # NetCDF/Exodus file ID
                         ex_entity_id  map_id,  # Elemental comm map ID
                         void_int* elem_ids,    # Element IDs
                         void_int* side_ids,    # Element side IDs
                         void_int* proc_ids,    # Processor IDs
                         int processor)         # This processor ID

    int ex_put_elem_cmap(int exoid,             # NetCDF/Exodus file ID
                         ex_entity_id  map_id,	# Elemental comm map ID
                         void_int* elem_ids,	# Vector of element IDs
                         void_int* side_ids,    # Vector of side IDs
                         void_int* proc_ids,    # Vector of processor IDs
                         int processor)	        # This processor ID

    #---
    # Deprecated functions
    #---

    # int ex_get_nodal_var(int exoid,
    #                      int time_step,
    #                      int nodal_var_index,
    #                      int64_t num_nodes,
    #                      void* nodal_var_vals)
    #
    # int ex_get_nodal_var_time(int exoid,
    #                           int nodal_var_index,
    #                           int64_t node_number,
    #                           int beg_time_step,
    #                           int end_time_step,
    #                           void* nodal_var_vals)
    #
    # # Use ex_get_concat_sets()
    # int ex_get_concat_node_sets(int exoid,
    #                             void_int* node_set_ids,
    #                             void_int* num_nodes_per_set,
    #                             void_int* num_df_per_set,
    #                             void_int* node_sets_node_index,
    #                             void_int* node_sets_df_index,
    #                             void_int* node_sets_node_list,
    #                             void* node_sets_dist_fact)
    #
    #
    # int ex_get_concat_side_sets(int exoid,
    #                             void_int* side_set_ids,
    #                             void_int* num_elem_per_set,
    #                             void_int* num_dist_per_set,
    #                             void_int* side_sets_elem_index,
    #                             void_int* side_sets_dist_index,
    #                             void_int* side_sets_elem_list,
    #                             void_int* side_sets_side_list,
    #                             void* side_sets_dist_fact)
    #
    # int ex_get_elem_attr(int exoid,
    #                      ex_entity_id elem_blk_id,
    #                      void* attrib)
    #
    # int ex_get_elem_attr_names(int exoid,
    #                            ex_entity_id elem_blk_id,
    #                            char** names)
    #
    # int ex_get_elem_blk_ids(int exoid,
    #                         void_int* ids)
    #
    # int ex_get_elem_block(int exoid,
    #                       ex_entity_id  elem_blk_id,
    #                       char* elem_type,
    #                       void_int* num_elem_this_blk,
    #                       void_int* num_nodes_per_elem,
    #                       void_int* num_attr)
    #
    # int ex_get_elem_conn(int exoid,
    #                      ex_entity_id elem_blk_id,
    #                      void_int* connect)
    #
    # int ex_get_elem_map(int exoid,
    #                     ex_entity_id map_id,
    #                     void_int* elem_map)
    #
    # int ex_get_elem_num_map(int exoid,
    #                         void_int* elem_map)
    #
    # int ex_get_elem_var(int exoid,
    #                     int time_step,
    #                     int elem_var_index,
    #                     ex_entity_id elem_blk_id,
    #                     int64_t num_elem_this_blk,
    #                     void* elem_var_vals)
    #
    # int ex_get_elem_var_tab(int exoid,
    #                         int num_elem_blk,
    #                         int num_elem_var,
    #                         int* elem_var_tab)
    #
    # int ex_get_elem_var_time(int exoid,
    #                          int elem_var_index,
    #                          int64_t elem_number,
    #                          int beg_time_step,
    #                          int end_time_step,
    #                          void* elem_var_vals)
    #
    # int ex_get_map(int exoid,
    #                void_int* elem_map)
    #
    # int ex_get_node_map(int exoid,
    #                     ex_entity_id map_id,
    #                     void_int* node_map)
    #
    # int ex_get_node_num_map(int exoid,
    #                         void_int* node_map)
    #
    # int ex_get_node_set_param(int exoid,
    #                           ex_entity_id  node_set_id,
    #                           void_int* num_nodes_in_set,
    #                           void_int* num_df_in_set)
    #
    # int ex_get_node_set(int exoid,
    #                     ex_entity_id node_set_id,
    #                     void_int* node_set_node_list)
    #
    # int ex_get_node_set_dist_fact (int exoid,
    #                                ex_entity_id node_set_id,
    #                                void* node_set_dist_fact)
    #
    # int ex_get_node_set_ids(int exoid,
    #                         void_int* ids)
    #
    # int ex_get_nset_var_tab(int exoid,
    #                         int num_nodesets,
    #                         int num_nset_var,
    #                         int* nset_var_tab)
    #
    # int ex_get_nset_var(int exoid,
    #                     int time_step,
    #                     int nset_var_index,
    #                     ex_entity_id nset_id,
    #                     int64_t num_node_this_nset,
    #                     void* nset_var_vals)
    #
    # int ex_get_one_elem_attr(int exoid,
    #                          ex_entity_id elem_blk_id,
    #                          int attrib_index,
    #                          void* attrib)
    #
    # int ex_get_side_set(int exoid,
    #                     ex_entity_id side_set_id,
    #                     void_int* side_set_elem_list,
    #                     void_int* side_set_side_list)
    #
    # int ex_get_side_set_dist_fact(int exoid,
    #                               ex_entity_id side_set_id,
    #                               void* side_set_dist_fact)
    #
    # int ex_get_side_set_ids(int exoid,
    #                         void_int* ss_ids)
    #
    # int ex_get_side_set_param(int exoid,
    #                           ex_entity_id  side_set_id,
    #                           void_int* num_side_in_set,
    #                           void_int* num_dist_fact_in_set)
    #
    # int ex_get_sset_var(int exoid,
    #                     int time_step,
    #                     int sset_var_index,
    #                     ex_entity_id sset_id,
    #                     int64_t num_side_this_sset,
    #                     void* sset_var_vals)
    #
    # int ex_get_sset_var_tab(int exoid,
    #                         int num_sidesets,
    #                         int num_sset_var,
    #                         int* sset_var_tab)
    #
    # int ex_get_var_names(int exoid,
    #                      const char* var_type,
    #                      int num_vars,
    #                      char* var_names[])
    #
    # int ex_get_var_name(int exoid,
    #                     const char* var_type,
    #                     int var_num,
    #                     char* var_name)
    #
    # int ex_get_var_param(int exoid,
    #                      const char* var_type,
    #                      int* num_vars)
    #
    # int ex_get_var_tab(int exoid,
    #                    const char* var_type,
    #                    int num_blk,
    #                    int num_var,
    #                    int* var_tab)
    #
    # int ex_put_concat_node_sets(int exoid,
    #                             void_int* node_set_ids,
    #                             void_int* num_nodes_per_set,
    #                             void_int* num_dist_per_set,
    #                             void_int* node_sets_node_index,
    #                             void_int* node_sets_df_index,
    #                             void_int* node_sets_node_list,
    #                             void* node_sets_dist_fact)
    #
    # int ex_put_concat_side_sets(int exoid,
    #                             void_int* side_set_ids,
    #                             void_int* num_elem_per_set,
    #                             void_int* num_dist_per_set,
    #                             void_int* side_sets_elem_index,
    #                             void_int* side_sets_dist_index,
    #                             void_int* side_sets_elem_list,
    #                             void_int* side_sets_side_list,
    #                             void* side_sets_dist_fact)
    #
    # int ex_put_concat_var_param(int exoid, int num_g, int num_n,
    #                             int num_e, int num_elem_blk, int* elem_var_tab)
    #
    # int ex_put_elem_attr_names(int exoid,
    #                            ex_entity_id elem_blk_id,
    #                            char* names[])
    #
    # int ex_put_elem_attr(int exoid,
    #                      ex_entity_id elem_blk_id,
    #                      const void* attrib)
    #
    # int ex_put_elem_block(int exoid,
    #                       ex_entity_id elem_blk_id,
    #                       const char* elem_type,
    #                       int64_t num_elem_this_blk,
    #                       int64_t num_nodes_per_elem,
    #                       int64_t num_attr)
    #
    # int ex_put_elem_conn(int exoid,
    #                      ex_entity_id elem_blk_id,
    #                      const void_int* connect)
    #
    # int ex_put_elem_map(int exoid,
    #                     ex_entity_id map_id,
    #                     const void_int* elem_map)
    #
    # int ex_put_elem_num_map(int exoid,
    #                         const void_int* elem_map)
    #
    # int ex_put_elem_var(int exoid,
    #                     int time_step,
    #                     int elem_var_index,
    #                     ex_entity_id elem_blk_id,
    #                     int64_t num_elem_this_blk,
    #                     const void* elem_var_vals)
    #
    # int ex_put_elem_var_tab(int exoid,
    #                         int num_elem_blk,
    #                         int num_elem_var,
    #                         int* elem_var_tab)
    #
    # int ex_put_glob_vars(int exoid,
    #                      int time_step,
    #                      int num_glob_vars,
    #                      const void* glob_var_vals)
    #
    # int ex_put_map(int exoid,
    #                const void_int* elem_map)
    #
    # int ex_put_node_map(int exoid,
    #                     ex_entity_id map_id,
    #                     const void_int* node_map)
    #
    # int ex_put_node_num_map(int exoid,
    #                         const void_int* node_map)
    #
    # int ex_put_node_set(int exoid,
    #                     ex_entity_id node_set_id,
    #                     const void_int* node_set_node_list)
    #
    # int ex_put_node_set_dist_fact (int exoid,
    #                                ex_entity_id node_set_id,
    #                                const void* node_set_dist_fact)
    #
    # int ex_put_node_set_param(int exoid,
    #                           ex_entity_id node_set_id,
    #                           int64_t num_nodes_in_set,
    #                           int64_t num_dist_in_set)
    #
    # int ex_put_nset_var(int exoid,
    #                     int time_step,
    #                     int nset_var_index,
    #                     ex_entity_id nset_id,
    #                     int64_t num_nodes_this_nset,
    #                     const void* nset_var_vals)
    #
    # int ex_put_nset_var_tab(int exoid,
    #                         int num_nset,
    #                         int num_nset_var,
    #                         int* nset_var_tab)
    #
    # int ex_put_one_elem_attr(int exoid,
    #                          ex_entity_id elem_blk_id,
    #                          int attrib_index,
    #                          const void* attrib)
    #
    # int ex_put_side_set(int exoid,
    #                     ex_entity_id side_set_id,
    #                     const void_int* side_set_elem_list,
    #                     const void_int* side_set_side_list)
    #
    # int ex_put_side_set_dist_fact(int exoid,
    #                               ex_entity_id side_set_id,
    #                               const void* side_set_dist_fact)
    #
    # int ex_put_side_set_param(int exoid,
    #                           ex_entity_id side_set_id,
    #                           int64_t num_side_in_set,
    #                           int64_t num_dist_fact_in_set)
    #
    # int ex_put_sset_var(int exoid,
    #                     int time_step,
    #                     int sset_var_index,
    #                     ex_entity_id sset_id,
    #                     int64_t num_faces_this_sset,
    #                     const void* sset_var_vals)
    #
    # int ex_put_sset_var_tab(int exoid,
    #                         int num_sset,
    #                         int num_sset_var,
    #                         int* sset_var_tab)
    #
    # int ex_put_var_name(int exoid,
    #                     const char* var_type,
    #                     int var_num,
    #                     const char* var_name)
    #
    # int ex_put_var_names(int exoid,
    #                      const char* var_type,
    #                      int num_vars,
    #                      char* var_names[])
    #
    # int ex_put_var_param(int exoid,
    #                      const char* var_type,
    #                      int num_vars)
    #
    # int ex_put_var_tab(int exoid,
    #                    const char* var_type,
    #                    int num_blk,
    #                    int num_var,
    #                    int* var_tab)

    #---
    # End of Deprecated functions and their replacements
    #---

    char* ex_name_of_object(ex_entity_type obj_type)
    
    ex_entity_type ex_var_type_to_ex_entity_type(char var_type)

    # Should be internal use only, but was in external include file for
    # nemesis and some codes are using the function
    int ex_get_idx(int neid,                # NetCDF/Exodus file ID
                   const char* ne_var_name, # Nemesis index variable name
                   int64_t *index,          # array of length 2 to hold results
                   int pos)                 # position of this proc/cmap in index