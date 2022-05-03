import React, { useEffect, useState } from 'react';
import ProjectService from 'shared/services/Project.service';
import ProjectList from 'components/projects/ProjectList';

function TeamProjectList({ team, forPart }) { // forPart can be "geneMapper" or "geneCruncher"
  const [projects, setProjects] = useState([]);
  const [loading, setLoading] = useState(true);

  useEffect(async () => {
    setLoading(true);
    setProjects(await ProjectService.getTeamProjects(team.id, forPart));
    setLoading(false);
  }, [team, forPart]);

  return <ProjectList isLoading={loading} projects={projects} />;
}

export default TeamProjectList;
