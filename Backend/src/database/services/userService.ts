import {IUser, userModel} from "../models/user";

export default class UserService {

    static async getAllUsers() {
        const allUsers = await userModel.find().exec();
        //const allUsers = await userModel.find();
        console.log("userService.getAllUsers : " + allUsers); // TODO
        return allUsers;
    }

    static async getUserById(user: any) {
        return userModel.findById( user._id );
    }

    static async getUserByEmail(email: string) {
        return userModel.findOne({email});
    }

    static async addUser(user: any): Promise<IUser> {
        console.log("userService.addUser : adding user " + user.toObject());
        let userAdded : (IUser | undefined) = undefined;
        userAdded = await userModel.create({
            firstName: user.firstName,
            lastName: user.lastName,
            email: user.email,
            password: user.saltHashedPassword,
            note: user.note
        });
        console.log("userService.addUser : " + userAdded.toObject()); // TODO
        return userAdded;
    }

    static async addUsers() {} // TODO

    static async updateUser() {} // TODO

    static async deleteUser() {} // TODO

    static async deleteUsers() {} // TODO
}
