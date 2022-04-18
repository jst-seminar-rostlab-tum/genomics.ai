import React, { useState, useEffect } from 'react';
import HeaderView from 'components/general/HeaderView';
import ProjectCard from 'components/projectOverview/ProjectCard';
import styles from './projectOverview.module.css';
import queryMyProjects from 'shared/services/mock/projects';

function ProjectOverview({ sidebarShown }) {
  const [projects, setProjects] = useState([]);
  useEffect(() => {
    queryMyProjects()
      .then((newProjects) => setProjects(newProjects))
      .catch((ignored) => { console.log(ignored); });
  }, [setProjects]);

  function onLeft(project) {
    setProjects(projects.filter((i) => i.id !== project.id));
  }

  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title="My Projects"
    >
      <div className={styles.content}>
        {projects.length === 0 ? 'No projects.' : ''}
        {projects.map((project) => (
          <div key={project.id}>
            <ProjectCard
              project={project}
              onLeft={(proj) => onLeft(proj)}
            />
            <div className={styles.cardSpacing} />
          </div>
        ))}
      </div>
    </HeaderView>
  );
}

export default ProjectOverview;
