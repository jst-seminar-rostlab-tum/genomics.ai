import axiosInstance from './axiosInstance';

const MODEL = 'users';

const UserService = {
  getUsers: async (params) => {
    const { data } = await axiosInstance.get(`/${MODEL}`, { params });
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
    return data;
  },
  getUserTeams: async (id) => {
    const { data } = await axiosInstance.get(`${MODEL}/${id}/teams`)
    return data;
  }
};

export default UserService;
