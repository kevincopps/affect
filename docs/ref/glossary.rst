========================================
Glossary
========================================

Exodus Database
---------------

Conceptual data model of the ExodusII database and :mod:`affect.exodus` module.

.. glossary::

   attributes
      Optional floating point numbers that can be assigned to each and every entry in a :class:`.Nodal`, :class:`.Set`,
      or :class:`.Block` entity. Every entry in an entity must have the same number of attributes, but the attribute
      values vary among the entries. Attributes are accessed through the member functions of an entity, for example,
      :meth:`.Set.num_attributes`, :meth:`.Set.attribute_names`, :meth:`.Set.attribute`, and :meth:`.Set.attributes`.

   block
      A association of entries with the same topology containing node connectivity information. For example, an
      element :class:`.Block` is an association of element entries and the nodes connected to each element.

   blocks
      A dictionary-like container of the :class:`.Block` instances of a certain :class:`.EntityType` in a
      :class:`.Database`. The dictionaries of different types of blocks are accessible through the following attributes:
      :attr:`.Database.edge_blocks`,
      :attr:`.Database.face_blocks`, or
      :attr:`.Database.element_blocks`.

   coordinates
      A special field associated with the :class:`.Nodal` entity storing the spatial coordinates of every node entry in
      the :class:`.Database`.

   database
      File storage for a mesh data model, a :class:`.Database` contains all the mesh entities and their corresponding
      entries, and the temporal :class:`.Field` variables.

   distribution factors
      Optional floating point values associated with every entry and every :class:`.Set` of a certain type, if they
      exist. Distribution factors are typically used in the simulation as a multiplier on an external load.
      Distribution factors are accessed through
      :meth:`.Sets.num_distribution_factors_all`,
      :meth:`.Set.num_distribution_factors`, and
      :meth:`.Set.distribution_factors`.

   entity
      An association of a subset of entries of a certain type (elements, faces, sides, edges, nodes). An entity is
      either the single :class:`.Global` or :class:`.Nodal` entity of the :class:`.Database`, or one of the possible
      multiple members of the :class:`.Blocks`, :class:`.Sets`, or :class:`.Maps` entities of the :class:`.Database`.

   entity ID
      An integer associated with each existing entity in the :class:`.Database`, the integer is unique to each entity
      of the same :class:`.EntityType`. The entity ID's are used as the keys used to access the dictionary-like
      containers :class:`.Blocks`, :class:`Sets`, and :class:`Maps`.

   entity Type
      One of the values of the enum :class:`.EntityType`, including 
      :data:`~.EntityType.NODAL`, 
      :data:`~.EntityType.NODE_SET`, 
      :data:`~.EntityType.EDGE_BLOCK`, 
      :data:`~.EntityType.EDGE_SET`, 
      :data:`~.EntityType.FACE_BLOCK`, 
      :data:`~.EntityType.FACE_SET`, 
      :data:`~.EntityType.ELEM_BLOCK`, 
      :data:`~.EntityType.ELEM_SET`, 
      :data:`~.EntityType.SIDE_SET`, 
      :data:`~.EntityType.ELEM_MAP`, 
      :data:`~.EntityType.NODE_MAP`, 
      :data:`~.EntityType.EDGE_MAP`, 
      :data:`~.EntityType.FACE_MAP`, 
      :data:`~.EntityType.GLOBAL`, and
      :data:`~.EntityType.COORDINATE`.

   entry
      Entries are the fundamental building blocks of the grid or mesh of a database. Entries refer to nodes, edges,
      faces, and elements of the :class:`.Database` mesh. Entries are not represented by their own Python objects,
      entry IDs, but they correspond to the first index of the :class:`.FieldArrays`.

   field
      A name for an array of values and the name of components associated with entries. The :class:`.Field` names are
      used to access the :class:`.FieldArray` values stored in the :class:`.Database`. Each of the named components of a
      :class:`.Field` with values in a :class:`.FieldArray` are a scalar :term:`variable` in the :class:`Database`.
      A :term:`field` is a grouping of ExodusII :term:`variable` by common name prefix;
      the suffix of the :term:`variable` name (whatever follows the last underscore '_' in the name) becomes a
      component name of the field. See also :term:`field array`.

   field array
      The actual scalar, vector and tensor values accessed in the :class:`.Database` by using a :class:`.Field` name and
      components. The :class:`.FieldArray` is a multidimensional array, with the first index corresponding to entries.
      It contains floating point values that vary in space (by :term:`entry` index) and time (:term:`time step`).
      Entities that may have field array values: global, nodal, blocks, and sets. For fields on blocks or sets, the
      field may or may not be active on all entities of that type; to find out use :meth:`.Block.is_field_active` or
      :meth:`.Set.is_field_active`. The values of the field array may be accessed on all entries at a single time step,
      for example see :meth:`.Nodal.field`; or on a range of entries at a time step, for example,
      :meth:`.Nodal.partial_field`; or on a single entry at all existing time steps, for example,
      :meth:`.Nodal.field_at_times`.

   global
      A :class:`.Global` is a single top level :class:`.Database` entity maintaining the spatial dimension,
      the number of time steps, the sums of all the entries of various types in the mesh (elements, faces, nodes)
      referenced in other :class:`.Database` entities. It is accessed from the attribute :attr:`.Database.global`.

   information data
      Info data is a list of optional supplementary text strings associated with a database. Typically this might be
      the input file from the simulation run that was executed to create the database results. Information data is
      accessed through :attr:`.Database.info`

   internal numbering
      The internal numbering of node entries is in the range [0, :meth:`.Global.num_nodes`]. The internal numbering of
      elements is by total subsequent entries in the :class:`.Block` in :meth:`.Database.blocks` (of
      type :data:`.EntityType.ELEMENT_BLOCK`) and these are in the range [0, :meth:`.Global.num_elements`].

   map
      A :class:`.Map` is a container of entries with new integers representing a number other than that of the
      default internal numbering for that type of entry.

   maps
      A dictionary-like container of the :class:`.Map` instances of a certain :class:`.EntityType` in a
      :class:`.Database`.
      The dictionaries of different types of maps are accessible through the following attributes:
      :attr:`.Database.element_maps`,
      :attr:`.Database.node_maps`,
      :attr:`.Database.edge_maps`, or
      :attr:`.Database.face_maps`.

   quality assurance records
      QA data are optional text strings in the :class:`.Database`, storing a history of application codes that modified
      the :class:`.Database` file, including the application name, description, date and time. Quality assurance data
      is accessed through :attr:`.Database.qa_records`

   nodal
      The single entity of a :class:`.Database` that stores nodal coordinates, nodal fields, and nodal attributes.
      The :class:`.Nodal` object is accessed from :attr:`.Database.nodal`.

   properties
      Optional named integer variables associated with every entity of a certain type in the database. The types of
      entities that may have properties are: :class:`.Block`, :class:`.Map`, and :class:`.Set` entities.
      Property names are accessed through the member function of the collection of entities, for example,
      :meth:`.Blocks.property_names`.
      Property values are accessed through the member functions of an entity, for example,
      :meth:`.Block.property`.

   set
      A :class:`.Set` entity is a container of a subset of the entries of a certain type (nodes, edges, faces, sides,
      elements) in the :class:`.Database`. There may be multiple sets of a certain type and they may intersect. Sets
      are usually used to apply boundary conditions to portions of the mesh, and sets may contain
      :term:`attributes`, :term:`properties` and :term:`distribution factors`, and multiple :term:`variable`.

   sets
      A dictionary-like container of the :class:`.Set` instances of a certain :class:`.EntityType` in a
      :class:`.Database`. The dictionaries of different types of sets are accessible through the following attributes:
      :attr:`.Database.node_sets`,
      :attr:`.Database.edge_sets`,
      :attr:`.Database.face_sets`, or
      :attr:`.Database.side_sets`. Entries of side sets are actually the pairing of an element and a local side number.

   variable
      Variables, in a :class:`.Database` are named scalar floating point arrays. The values of variables vary in time
      and are associated with entries in the database. A single variable is one component of a
      more useful multi-dimensional :class:`FieldArray`, there is often no need to refer to variables
      separately from a :class:`FieldArray`. The suffix of a name of a Exodus variable is also the name of a
      :class:`Field` component. The underlying scalar variable values making up :term:`field array` may be accessed
      in the database in a similar way to their :class:`FieldArray` counterpart.

   time step
      The discrete values of time at which the values of fields (variables) are stored in the database. The values of
      time steps are accessible through the attribute :attr:`Database.globals.num_times` and
      :attr:`Database.globals.times`.
