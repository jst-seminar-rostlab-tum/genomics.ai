type ImageUploadResult = ImageUploadResultOk|ImageUploadResultError

interface ImageUploadResultOk{
    success:true
    objectUrl:string
}

interface ImageUploadResultError {
    success: false,
    status: number,
    message: string,
    error?: any
}
