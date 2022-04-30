import axiosInstance from "./axiosInstance"

const MODEL = "jobs"

const JobService = {
    getJobs: async () => {
        const { data } = await axiosInstance.get(`/${MODEL}`)
        return data
    },
}

export default JobService