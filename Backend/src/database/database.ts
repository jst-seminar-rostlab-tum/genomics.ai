import mongoose from "mongoose";

export class Database {
  async connect() {
    const DATABASE_URI = process.env.DATABASE_URI;
    await mongoose.connect(DATABASE_URI!).then(
      () => console.log("Connected to database"),
      (err) => console.error(err)
    );
    return;
  }
}
