import axiosInstance from './axiosInstance';

const CONTACT_US = 'contactus';

const ContactUsService = {
  postContactForm: async ( data ) => {
    await axiosInstance.post(`/${CONTACT_US}`, data)
    .then((response)=>console.log(response))
    .catch((error)=>console.log(error))
  }
};

export default ContactUsService;
