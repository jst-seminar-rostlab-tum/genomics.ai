import React, { useEffect, useState } from 'react';
import ProjectService from 'shared/services/Project.service';
import ProjectBarCard from 'components/GeneMapper/projectBarCard';
import AtlasService from 'shared/services/Atlas.service';
import ModelService from 'shared/services/Model.service';
import TeamService from 'shared/services/Team.service';
import { useSubmissionProgress } from 'shared/context/submissionProgressContext';
import { MULTIPART_UPLOAD_STATUS, PROJECTS_UPDATE_INTERVAL, statusIsError } from 'shared/utils/common/constants';

function TeamProjectList({ team }) {
  const [projects, setProjects] = useState([]);
  const [atlases, setAtlases] = useState([]);
  const [models, setModels] = useState([]);
  const [userTeams, setUserTeams] = useState([]);
  const [submissionProgress, setSubmissionProgress] = useSubmissionProgress();

  console.log(projects, team._id);

  useEffect(() => {
    AtlasService.getAtlases().then((data) => setAtlases(data));
    ModelService.getModels().then((data) => setModels(data));
    TeamService.getMyTeams().then((teams) => setUserTeams(teams));
  }, []);

  useEffect(() => {
    ProjectService.getTeamProjects(team._id).then((data) => { console.log(data); setProjects(data); });
    if (submissionProgress.status === MULTIPART_UPLOAD_STATUS.COMPLETE
      || submissionProgress.status === MULTIPART_UPLOAD_STATUS.CANCELING
      || statusIsError(submissionProgress.status)) {
      setSubmissionProgress({
        status: MULTIPART_UPLOAD_STATUS.IDLE,
        uploadId: '',
        chunks: 0,
        uploaded: 0,
        remaining: [],
        uploadedParts: [],
      });
    }
  }, [submissionProgress.status]);

  return projects
    .map((project) => {
      console.log("im here!");
      return (
        <ProjectBarCard
          key={project._id}
          project={project}
          atlas={atlases.find((atlas) => String(atlas._id) === String(project.atlasId))}
          model={models.find((model) => String(model._id) === String(project.modelId))}
          userTeams={userTeams}
          addProjectToTeam={() => true}
          handleDelete={() => true}
          submissionProgress={submissionProgress.uploadId === project.uploadId
            ? submissionProgress : null}
          setSubmissionProgress={submissionProgress.uploadId === project.uploadId
            ? setSubmissionProgress : () => { }}
        />

      );
    })
}

export default TeamProjectList;
