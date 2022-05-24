import React, { useEffect, useState } from 'react';
import { CircularProgress } from '@mui/material';
import ProjectService from 'shared/services/Project.service';
import ProjectBarCard from 'components/GeneMapper/projectBarCard';
import AtlasService from 'shared/services/Atlas.service';
import ModelService from 'shared/services/Model.service';
import TeamService from 'shared/services/Team.service';

function TeamProjectList({ teamId }) {
  const [projects, setProjects] = useState([]);
  const [atlases, setAtlases] = useState([]);
  const [models, setModels] = useState([]);
  const [userTeams, setUserTeams] = useState([]);
  const [isLoading, setIsLoading] = useState(true);

  useEffect(() => {
    setIsLoading(true);
    AtlasService.getAtlases().then((data) => setAtlases(data));
    ModelService.getModels().then((data) => setModels(data));
    TeamService.getMyTeams().then((teams) => setUserTeams(teams));
    ProjectService.getTeamProjects(teamId).then((data) => setProjects(data));
    setIsLoading(false);
  }, []);

  if (isLoading) {
    return <CircularProgress />;
  }

  if (projects.length === 0) {
    return (<span>No projects.</span>);
  }

  return projects
    .map((project) => (
      <ProjectBarCard
        key={project._id}
        project={project}
        atlas={atlases.find((atlas) => String(atlas._id) === String(project.atlasId))}
        model={models.find((model) => String(model._id) === String(project.modelId))}
        userTeams={userTeams}
        handleDelete={() => true}
      />
    ));
}

export default TeamProjectList;
