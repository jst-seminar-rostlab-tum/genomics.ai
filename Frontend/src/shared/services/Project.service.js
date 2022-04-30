import axiosInstance from "./axiosInstance"

// To be consistent backend team should rename this to projects
const MODEL = "jobs"

const ProjectService = {
    getProjects: async () => {
        const { data } = await axiosInstance.get(`/${MODEL}`)
        return data
    },
}

export default ProjectService