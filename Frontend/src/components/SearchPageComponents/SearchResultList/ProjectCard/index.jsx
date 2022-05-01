import React from 'react';

import { Chip } from '@mui/material';
import SearchCard from '../SearchCard';
import LabeledLink from '../LabeledLink';

// Card to display search result for a single project
function ProjectCard({ item: project }) {
  return (
    <SearchCard
      title={project.name}
      link={`/sequencer/projects/${project.id}`}
      primary={<Chip label={project.type} color="primary" size="small" />}
      tertiary={(
        <>
          {project.institution && (
            <LabeledLink
              label="Institution"
              content={project.institution.name}
              to={`/sequencer/institutions/${project.institution.id}`}
            />
          )}
          {project.team && (
            <LabeledLink
              label="Team"
              content={project.team.name}
              to={`/sequencer/teams/${project.team.id}`}
            />
          )}
        </>
      )}
    />
  );
}

export default ProjectCard;
