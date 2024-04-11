from icecube import dataclasses


def combine_i3_particle(
    frame,
    output_name,
    pos_x_name=None,
    pos_y_name=None,
    pos_z_name=None,
    dir_name=None,
    time_name=None,
    energy_name=None,
    shape="Cascade",
):
    """Add an I3Particle to the frame with values from other particles.

    This function can be helpful when creating seeds for reconstruction
    methods based on previous reconstruction results.

    Parameters
    ----------
    frame : I3Frame
        The current physics frame
    output_name : str
        The frame key to which the combined I3Particle will be written to.
    pos_x_name : None, optional
        The name of the I3Particle from which to take the vertex-x position.
    pos_y_name : None, optional
        The name of the I3Particle from which to take the vertex-y position.
    pos_z_name : None, optional
        The name of the I3Particle from which to take the vertex-z position.
    dir_name : None, optional
        The name of the I3Particle from which to take the direction.
    time_name : None, optional
        The name of the I3Particle from which to take the time of the vertex.
    energy_name : None, optional
        The name of the I3Particle from which to take the energy.
    shape : str, optional
        The dataclasses.I3Particle.ParticleShape to assign.
    """
    particle = dataclasses.I3Particle()
    if pos_x_name is not None:
        particle.pos.x = frame[pos_x_name].pos.x
    if pos_y_name is not None:
        particle.pos.y = frame[pos_y_name].pos.y
    if pos_z_name is not None:
        particle.pos.z = frame[pos_z_name].pos.z
    if dir_name is not None:
        particle.dir = frame[dir_name].dir
    if time_name is not None:
        particle.time = frame[time_name].time
    if energy_name is not None:
        particle.energy = frame[energy_name].energy

    particle.fit_status = dataclasses.I3Particle.FitStatus.OK
    particle.shape = getattr(dataclasses.I3Particle.ParticleShape, shape)

    frame[output_name] = particle


def combine_i3_particle_from_dict(
    frame,
    output_name,
    pos_x_name=None,
    pos_y_name=None,
    pos_z_name=None,
    dir_names=None,
    time_name=None,
    energy_name=None,
    shape="Cascade",
):
    """Add an I3Particle to the frame with values from other dict keys,
    that are not I3Particles. The specific key of the dict must be
    separated by a dot. For example: "RecoPosition.center_pos_x".

    This function can be helpful when creating seeds for reconstruction
    methods based on previous reconstruction results.

    Parameters
    ----------
    frame : I3Frame
        The current physics frame
    output_name : str
        The frame key to which the combined I3Particle will be written to.
    pos_x_name : None, optional
        The name of the dict from which to take the vertex-x position.
    pos_y_name : None, optional
        The name of the dict from which to take the vertex-y position.
    pos_z_name : None, optional
        The name of the dict from which to take the vertex-z position.
    dir_names : list, optional
        The name of the dict from which to take the direction. It is possible
        to provide 2 or 3 keys, depending on the dimension of the direction.
        2 keys: (zenith, azimuth),
        3 keys: (x, y, z).
    time_name : None, optional
        The name of the dict from which to take the time of the vertex.
    energy_name : None, optional
        The name of the dict from which to take the energy.
    shape : str, optional
        The dataclasses.I3Particle.ParticleShape to assign.
    """

    particle = dataclasses.I3Particle()
    if pos_x_name is not None:
        particle.pos.x = frame[".".join(pos_x_name.split(".")[:-1])][
            pos_x_name.split(".")[-1]
        ]
    if pos_y_name is not None:
        particle.pos.y = frame[".".join(pos_y_name.split(".")[:-1])][
            pos_y_name.split(".")[-1]
        ]
    if pos_z_name is not None:
        particle.pos.z = frame[".".join(pos_z_name.split(".")[:-1])][
            pos_z_name.split(".")[-1]
        ]
    if dir_names is not None:
        if len(dir_names) == 2:
            particle.dir = dataclasses.I3Direction(
                frame[".".join(dir_names[0].split(".")[:-1])][
                    dir_names[0].split(".")[-1]
                ],
                frame[".".join(dir_names[1].split(".")[:-1])][
                    dir_names[1].split(".")[-1]
                ],
            )
        elif len(dir_names) == 3:
            particle.dir = dataclasses.I3Direction(
                frame[".".join(dir_names[0].split(".")[:-1])][
                    dir_names[0].split(".")[-1]
                ],
                frame[".".join(dir_names[1].split(".")[:-1])][
                    dir_names[1].split(".")[-1]
                ],
                frame[".".join(dir_names[2].split(".")[:-1])][
                    dir_names[2].split(".")[-1]
                ],
            )
        else:
            raise ValueError("dir_names must be of length 2 or 3", dir_names)
    if time_name is not None:
        particle.time = frame[".".join(time_name.split(".")[:-1])][
            time_name.split(".")[-1]
        ]
    if energy_name is not None:
        particle.energy = frame[".".join(energy_name.split(".")[:-1])][
            energy_name.split(".")[-1]
        ]

    particle.fit_status = dataclasses.I3Particle.FitStatus.OK
    particle.shape = getattr(dataclasses.I3Particle.ParticleShape, shape)

    frame[output_name] = particle
