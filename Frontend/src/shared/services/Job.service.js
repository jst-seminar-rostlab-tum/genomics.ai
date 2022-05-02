import axiosInstance from './axiosInstance';
import MockJobService from './mock/Job.service';

const MOCK_JOBS = true;

const JobService = MOCK_JOBS ? MockJobService : {
  getJobs: async () => {
    const { data } = await axiosInstance.get('/jobs');
    return data;
  },
};

export default JobService;
