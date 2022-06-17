import React from 'react';
import LoadingList from 'components/general/LoadingList';
import ProjectCard from 'components/projects/ProjectCard';

function ProjectList({ isLoading, projects }) {
  return (
    <LoadingList
      isLoading={isLoading}
      elements={projects}
      cardBuilder={(project) => <ProjectCard project={project} />}
      noElementsMessage="No projects."
    />
  );
}

export default ProjectList;
