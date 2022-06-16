import axiosInstance from './axiosInstance';

const MODEL = 'users';

function enhanceInstitution(institution) {
  return { ...institution, id: institution._id };
}

function enhanceTeam(team) {
  return { ...team, id: team._id, name: team.title };
}


const UserService = {
  getUsers: async (params) => {
    const preparedParams = { ...params };
    if (!preparedParams.sortBy || preparedParams.sortBy === 'name') {
      preparedParams.sortBy = 'lastName';
    }
    const { data } = await axiosInstance.get(`/${MODEL}`, { params: preparedParams });
    return data;
  },

  getOwnTeams: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}/ownteams`);
    return data;
  },

  getUser: async (id) => {
    const { data } = await axiosInstance.get(`${MODEL}/${id}`)
    return data;
  },
  getUserInstitutions: async (id) => {
    const { data } = await axiosInstance.get(`${MODEL}/${id}/institutions`)
    return data.map(enhanceInstitution);;
  },
  getUserTeams: async (id) => {
    const { data } = await axiosInstance.get(`${MODEL}/${id}/teams`)
    return data.map(enhanceTeam);
  }
};

export default UserService;
