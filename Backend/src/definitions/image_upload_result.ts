export type ImageUploadResult = ImageUploadResultOk|ImageUploadResultError

export interface ImageUploadResultOk{
    success:true
    objectUrl:string
}

export interface ImageUploadResultError {
    success: false,
    status: number,
    message: string,
    error?: any
}
