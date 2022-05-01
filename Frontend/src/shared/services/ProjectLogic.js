import axios from 'axios';
import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
import { getAuthAndJsonHeader, getAuthHeader } from 'shared/utils/common/utils';
import { startOrContinueUpload } from './UploadLogic';

function api(endpoint) {
  return `${BACKEND_ADDRESS}${endpoint.startsWith('/') ? '' : '/'}${endpoint}`;
}

async function get(endpoint) {
  return axios.get(api(endpoint), { headers: getAuthHeader() })
    .then((response) => response.data);
}

export default async function getProjects() {
  return get('/projects');
}

export async function getProject(id) {
  return get(`/project/${id}`);
}

export async function startOrContinueProjectUpload(
  selectedFile,
  submissionProgress,
  setSubmissionProgress,
  projectData,
) {
  return startOrContinueUpload(selectedFile,
    submissionProgress,
    setSubmissionProgress,
    projectData);
}

export async function addProjectToTeam(teamId, projectId) {
  return axios.put(api(`/teams/${teamId}/add_project`), { projectId },
    { headers: getAuthAndJsonHeader() });
}

export async function getOwnTeams() {
  return get('/users/ownteams');
}

export async function getModel(id) {
  return get(`/model/${id}`);
}

export async function getAtlas(id) {
  return get(`/atlas/${id}`);
}

export async function getModels() {
  return get('/models');
}

export async function getAtlases() {
  return get('/atlases');
}
