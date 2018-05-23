Calculation on all timesteps by block
=====================================

Here is a relatively larger example where you want to perform a operation on
every time step, on all element blocks, involving all the nodal fields in the
database. Suppose your calculation is performed in the function
``my_block_calculation``. Just for the sake of this example, we suppose that
the result of the block calculation is a scalar, and that these are summed into
a global result, which is further summed over time steps.

Outline
-------

The general idea of the procedure is,

* open the database
    * get the node coordinates
    * for each block
        * get the element-to-node connectivity and keep in memory
    * for each time step
        * get global field values
        * for each block
            * gather the field values to local values on the block
            * perform your calculation

Here since we want to operate on all the nodal fields block-by-block, for each
time step first read in all the nodal field arrays. Following that, as we loop
over blocks, we create an indexed local field array for a block.

Our strategy uses extra memory for local copies of the field values as we
operate on a block, but makes the calculation efficient in terms of speed. The
memory for the arrays storing the local copies of field values for one block
can be garbage collected as you proceed to work on the next block.

Keeping local connectivities
----------------------------

The local connectivities (element-to-node ID connectivity for a block) are
kept in compressed form in memory until the point at which they are needed.
This is handled transparently by the
:class:`affect.exodus.LocalConnectivity` and using the ``compress=True``
option passed to :meth:`affect.exodus.Blocks.connectivity_local_all`,
emphasized below. After the connectivity is used in uncompressed form, the
uncompressed copy can be garbage collected.

Copying fields from global to local
-----------------------------------

The local indexed copying of the nodal coordinates and other nodal field arrays
are created by using the :meth:`numpy.ndarray.take`, emphasized below.

.. code-block:: python
    :caption: block_by_block_calculation_example.py
    :name: block_by_block_calculation_example
    :linenos:
    :emphasize-lines: 14,19-20,46

    with exodus.DatabaseFile('/tmp/myExodusFile.exo') as e:

        num_times = e.globals.num_times()  # read number of time steps

        nodal = e.nodal                    # get the nodal object
        fields = exodus.Fields(nodal)      # ordered dictionary of field info
        coordinates = nodal.coordinates()  # read array for all nodes

        # Read the block local connectivities and store in a dictionary.
        # And the same for node coordinates, since they don't change with
        # time step.
        local_connectivities = OrderedDict()  # maintain block order
        local_coordinates = dict()
        local_iterator = e.element_blocks.connectivity_local_all(compress=True)
        for block_id, block, local in local_iterator:
            local_connectivities[block_id] = local
            # We use the take function to select the
            # global nodes to copy into our local array
            local_coordinates[block_id] = coordinates.take(local.global_nodes,
                                                            axis=0)

        all_times_result = 0.0

        for time_step in range(num_times):  # loop over time steps

            # Read the value of all nodal field arrays on global nodes at this
            # time step. Here you may decide to select only the subset of
            # fields you need for your calculation by name.
            global_arrays = dict()  # to hold the field arrays
            for name, field in fields:
                global_arrays[name] = nodal.field(field, time_step)  # read

            all_blocks_result = 0.0

            for block_id, local in local_connectivities:  # each block

                if local is None:
                    continue  # skip blocks without any nodes

                # Copy relevant node field values from global to
                # local arrays for this block.
                local_fields = dict()
                # Add the local coordinates first, which we already have.
                local_fields['coordinate'] = local_coordinates[block_id]
                for name, array in global_arrays:
                    local_fields[name] = array.take(local.global_nodes, axis=0)

                # Perform our calculation on local_fields arrays on this block.
                block_result = my_block_calculation(block_id,
                    local.local_nodes, local_fields)

                all_blocks_result += block_result

            all_times_result += all_blocks_result

So all the arrays passed to ``my_block_calculation`` are in local form
specific to a single block at a time.
